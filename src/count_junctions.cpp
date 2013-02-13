#include "count_junctions.h"
//#include <sstream>
#include <string>
#include <map>
#include <set>

using namespace std;

/*! @typedef
  @abstract Structure to provide the data to the bam_fetch() function
  @field junctionsU   map<string,int> with junction counts (non-allelic or allelic Unknown)
                      key of the form "chromosome:first_intronic_base:last_intronic_base:strand"
  @field junctionsR   map<string,int> with junction counts (allelic Reference)
  @field junctionsA   map<string,int> with junction counts (allelic Alternative)
  @field tname        target name string
  @field allelic      allelic true(1) or false(0)
 */
typedef struct {
    map<string,int> junctionsU;
    map<string,int> junctionsR;
    map<string,int> junctionsA;
    char* tname;
    int allelic;
} fetch_param;


/*! @function
  @abstract  callback for bam_fetch(); sums alignments if correct strand and shifted biological start/end position overlap fetch region
  @param  hit   the alignment
  @param  data  user provided data
  @return       0 if successful
 */
static int _addJunction(const bam1_t *hit, void *data){
    static uint8_t *xv_ptr = 0;
    fetch_param *jinfo = (fetch_param*)data;

    //ostringstream ss;
    //static string jid; // junction identifier, of the form "chromosome:first_intronic_base:last_intronic_base:strand"
    static char strbuffer[1024];
    static map<string,int>::iterator it;

    uint32_t *cigar = bam1_cigar(hit);
    int i = 0; // current position of cigar operations
    int l; // length of current cigar operation
    int x; // leftmost coordinate of the current cigar operation in reference string
    int y; // leftmost coordinate of the current cigar operation in read
    int op; // current cigar operation type

    // get pointer to junctions counter
    static map<string,int> *junctions = NULL;
    if(jinfo->allelic==0) {
	junctions = &(jinfo->junctionsU);

    } else {
	// get XV tag
	xv_ptr = bam_aux_get(hit,"XV");
	if(xv_ptr == 0)
	    Rf_error("XV tag missing but needed for allele-specific counting");
	switch(bam_aux2A(xv_ptr)){
	case 'U':
	    junctions = &(jinfo->junctionsU);
	    break;
	case 'R':
	    junctions = &(jinfo->junctionsR);
	    break;
	case 'A':
	    junctions = &(jinfo->junctionsA);
	    break;
	default:
	    Rf_error("'%c' is not a valid XV tag value; should be one of 'U','R' or 'A'", bam_aux2A(xv_ptr));
	}
    }
    
    // loop over cigar operations
    x = hit->core.pos;
    for (i = y = 0; i < hit->core.n_cigar; ++i) {

	l = cigar[i]>>4, op = cigar[i]&0xf;

	if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
	    // aligned bases -> increase read and reference coordinates
	    x += l; y += l;

	} else if (op == BAM_CREF_SKIP) { // N skipped region spliced alignment
	    // read skips reference region -> store junction and increase reference coordinate
	    sprintf(strbuffer, "%s:%i:%i:%c", jinfo->tname, x+1, x+l, hit->core.flag & BAM_FREVERSE ? '-' : '+');
	    //ss << jinfo->tname << ":" << (x+1) << ":" << (x+l) << ":" << (hit->core.flag & BAM_FREVERSE ? "-" : "+");
	    //jid = ss.str();
	    string jid(strbuffer);
	    it = junctions->find(jid);
	    if(it == junctions->end()) {
		junctions->insert(pair<string,int>(jid, 1));
	    } else {
		it->second++;
	    }
	    x += l;

	} else if (op == BAM_CDEL) {
	    // read skips reference region -> store junction and increase reference coordinate
	    sprintf(strbuffer, "%s:%i:%i:%c", jinfo->tname, x+1, x+l, hit->core.flag & BAM_FREVERSE ? '-' : '+');
	    //ss << jinfo->tname << ":" << (x+1) << ":" << (x+l) << ":" << (hit->core.flag & BAM_FREVERSE ? "-" : "+");
	    //jid = ss.str();
	    string jid(strbuffer);
	    it = junctions->find(jid);
	    if(it == junctions->end()) {
		junctions->insert(pair<string,int>(jid, 1));
	    } else {
		it->second++;
	    }
	    x += l;

	} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
	    // read has an insertion to the reference, or is soft-clipped  -> increase read coordinate
	    y += l;
	}
    }

    return 0;
}


/*! @function
  @abstract  enumerates and counts all introns (alignment-insertions) observed in a bam file
  @param   bamfile   Name of the bamfile
  @param   tid       integer vector of target identifiers
  @param   start     integer vector of target region start coordinates
  @param   end       integer vector of target region end coordinates
  @param   allelic   logical(1) to indicate allelic/non-allelic counting
  @return            allelic==FALSE: named vector of the junctions counts (names of the form "chromosome:first_intronic_base:last_intronic_base:strand")
                     allelic==TRUE: list of length 4 (names, R, U and A junction counts)
 */
SEXP count_junctions(SEXP bamfile, SEXP tid, SEXP start, SEXP end, SEXP allelic) {
    // check parameters
    if(!Rf_isString(bamfile) || Rf_length(bamfile) != 1)
        Rf_error("'bamfile' must be of type character(1)");
    if(!Rf_isInteger(tid))
        Rf_error("'tid' must be of type integer");
    if(!Rf_isInteger(start))
        Rf_error("'start' must be of type integer");
    if(!Rf_isInteger(end))
        Rf_error("'end' must be of type integer");
    if(!Rf_isLogical(allelic) || Rf_length(allelic) != 1)
        Rf_error("'allelic' must be of type logical(1)");

    // open bam file
    samfile_t *fin = 0;
    fin = samopen(Rf_translateChar(STRING_ELT(bamfile, 0)), "rb", NULL);
    if (fin == 0)
        Rf_error("failed to open BAM file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (fin->header == 0 || fin->header->n_targets == 0) {
        samclose(fin);
        Rf_error("BAM header missing or empty of file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }
    // open bam index
    bam_index_t *idx = 0; 
    idx = bam_index_load(Rf_translateChar(STRING_ELT(bamfile, 0)));
    if (idx == 0){
        samclose(fin);
        Rf_error("failed to open BAM index file: '%s'", Rf_translateChar(STRING_ELT(bamfile, 0)));
    }

    // initialise junction map
    fetch_param jinfo;
    if(Rf_asLogical(allelic))
	jinfo.allelic = 1;
    else
	jinfo.allelic = 0;

    // select bam_fetch callback function
    bam_fetch_f fetch_func = _addJunction;
    
    // loop over query regions
    int num_regions = Rf_length(tid);
    for(int i = 0; i < num_regions; i++){
	// set tname
	jinfo.tname = fin->header->target_name[INTEGER(tid)[i]];

        // process alignments that overlap region
        bam_fetch(fin->x.bam, idx, INTEGER(tid)[i], 
                  INTEGER(start)[i], // 0-based inclusive start
                  INTEGER(end)[i],   // 0-based exclusive end
                  &jinfo, fetch_func);
    }

    // clean up
    samclose(fin);
    bam_index_destroy(idx);

    SEXP count,  attrib;
    int n, i;
    if(jinfo.allelic == 1) {
	// store results in list SEXP
	//    get complete set of unique junctions
	set<string> uniqueJunctions;
	for(map<string,int>::iterator it=jinfo.junctionsU.begin(); it!=jinfo.junctionsU.end(); ++it)
	    uniqueJunctions.insert(it->first);
	for(map<string,int>::iterator it=jinfo.junctionsR.begin(); it!=jinfo.junctionsR.end(); ++it)
	    uniqueJunctions.insert(it->first);
	for(map<string,int>::iterator it=jinfo.junctionsA.begin(); it!=jinfo.junctionsA.end(); ++it)
	    uniqueJunctions.insert(it->first);
	n = uniqueJunctions.size();
	SEXP junctionNames, countR, countU, countA;
	PROTECT(junctionNames = allocVector(STRSXP, n));
	PROTECT(countR = allocVector(INTSXP, n));
	PROTECT(countU = allocVector(INTSXP, n));
	PROTECT(countA = allocVector(INTSXP, n));

	//    allocate and fill in count vectors
	i=0;
	map<string,int>::iterator cit;
	for(set<string>::iterator it=uniqueJunctions.begin(); it!=uniqueJunctions.end(); ++it) {
	    // id
	    SET_STRING_ELT(junctionNames, i, mkChar((*it).c_str()));
	    // R
	    cit = jinfo.junctionsR.find(*it);
	    INTEGER(countR)[i] = (cit == jinfo.junctionsR.end()) ? 0 : cit->second;
	    // U
	    cit = jinfo.junctionsU.find(*it);
	    INTEGER(countU)[i] = (cit == jinfo.junctionsU.end()) ? 0 : cit->second;
	    // A
	    cit = jinfo.junctionsA.find(*it);
	    INTEGER(countA)[i] = (cit == jinfo.junctionsA.end()) ? 0 : cit->second;

	    i++;
	}

	//    build list
	PROTECT(count = allocVector(VECSXP, 4));
	PROTECT(attrib = allocVector(STRSXP, 4));

	SET_STRING_ELT(attrib, 0, mkChar("id"));
	SET_STRING_ELT(attrib, 1, mkChar("R"));
	SET_STRING_ELT(attrib, 2, mkChar("U"));
	SET_STRING_ELT(attrib, 3, mkChar("A"));
	SET_VECTOR_ELT(count, 0, junctionNames);
	SET_VECTOR_ELT(count, 1, countR);
	SET_VECTOR_ELT(count, 2, countU);
	SET_VECTOR_ELT(count, 3, countA);
	setAttrib(count, R_NamesSymbol, attrib);
	UNPROTECT(6);

    } else {
	// store results in named vector SEXP
	n = jinfo.junctionsU.size();
	i=0;
	PROTECT(count = allocVector(INTSXP, n));
	PROTECT(attrib = allocVector(STRSXP, n));
	for(map<string,int>::iterator it=jinfo.junctionsU.begin(); it!=jinfo.junctionsU.end(); ++it) {
	    INTEGER(count)[i] = it->second;
	    SET_STRING_ELT(attrib, i, mkChar(it->first.c_str()));
	    i++;
	}
	setAttrib(count, R_NamesSymbol, attrib);
	UNPROTECT(2);
    }

    return count;
}
