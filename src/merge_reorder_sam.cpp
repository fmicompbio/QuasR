#include "merge_reorder_sam.h"
//#include "merge_reorder_sam_standalone.h"
#include <iostream>
#include <sstream>
#include <map>

using namespace std;
#define MAX_NM 10000 // nm tag value if read is not mapped 

class idLine { // stores a single alignment with integer identifier, flag and boolean isMapped 
public:
    int id;        // integer id prefix
    int bisQueue;  // which queue did the alignment come from (0..3)
    bool isMapped; // read(-pair) mapped?
    string line;   // SAM line (first read)
    string line2;  // SAM line (second read)
    idLine() {id=-1; bisQueue=-1; isMapped=false; line=""; line2=""; }
    idLine(const int &newid, const int &newbisQueue, const bool newIsMapped, const string &newline, const string &newline2)
    {id=newid; bisQueue=newbisQueue; isMapped=newIsMapped; line=newline; line2=newline2; }
    idLine(const int &newid, const bool newIsMapped, const string &newline, const string &newline2)
    {id=newid; bisQueue=-1; isMapped=newIsMapped; line=newline; line2=newline2; }
    bool operator() (const idLine& lhs, const idLine&rhs) const {return (lhs.id>rhs.id);}
    void print() { cerr << "  " << id << ":" << line; if(line2 != "") cerr << "; " << line2; cerr << endl; }
};

void _reverse_complement(string&);
void _replace_sequence(string&, bool);
void _remove_MD_tag(string&);
void _fix_FLAGs_and_sequences(idLine&);
int flush_bisulfite(int, ofstream&, map<int, string>&, vector<idLine>&, bool);  // same as flush_simple, bisulfite-version
int flush_allele(int, ofstream&, map<int, string>&, idLine&, char);     // same as flush_simple, allele-specific-version
int _make_unmapped_alignment(int, idLine&, map<int, string>&, bool, bool);
int _get_nm_tag(const idLine&);

class samFile { // handles a sam file
    static int nTotal; // number of samFile instances created
    static int nEof;   // number of samFile instances that reached fh.eof()

    const char *fname; // file name
    ifstream fh;       // input file stream

    string readbuffer; // current alignment line
    string readbuffer2;// current alignment line2 (paired reads)
    int readid;        // current alignment identifier
    bool readIsMapped; // current alignment is mapped
    bool readIsPaired; // current alignment is paired

    priority_queue<idLine, vector<idLine>, idLine> queue; // stores alignment until .flush()

    int getNextAln(); // read next (pair of) alignment, extract readid and flag
public:
    samFile(const char*);
    ~samFile();
    int advance(int);
    int flush_simple(int, ofstream&, map<int, string>&);           // output alignments to ofstream, store unmapped in map<>
    bool isEmpty() { return queue.empty(); }

    static bool allEof() { return (nTotal==nEof); }
    static int flush_unmapped(int, ofstream&, map<int, string>&, int);
    int get_nm_tag(int);    // get the sum of the nm-tags from the top of the queue
    int get_alignments_bisulfite(int, int, vector<idLine>&, map<int, string>&, bool, bool);
    int get_alignments_allele(int, vector<idLine>&, map<int, string>&);
};

int samFile::nTotal=0;
int samFile::nEof=0;

// constructor
samFile::samFile (const char* myfname) {
    fname = myfname;

    // open file
    fh.open(fname, ifstream::in | ifstream::binary);
    if(! fh.good()) {
	Rf_error("error opening file '%s'\n",fname);
    } else {
	// skip header
	while(fh.peek()=='@' && fh.good())
	    fh.ignore(INT_MAX, '\n');
    }

    nTotal++;
}

// destructor
samFile::~samFile () {
    if(fh.is_open())
        fh.close();
}


// read next (pair of) alignment, extract readid and flag
int samFile::getNextAln() {
    static size_t start_pos, end_pos;
    static int readflag, readid2, readflag2;
    static bool readIsMapped2;

    // read line
    getline (fh, readbuffer, '\n');
    if(fh.eof()) {
	nEof++;
	return 1;
    } else if(!fh.good()) {
	Rf_error("error reading from %s\n", fname);
    }

    // remove \r if exists (for windows)
    if(readbuffer[readbuffer.size()-1] == '\r')
        readbuffer.erase(readbuffer.size()-1, 1);

    // extract id
    end_pos = readbuffer.find('_');
    if(end_pos != string::npos)
	readid = atoi(readbuffer.substr(0, end_pos).c_str());
    else
	Rf_error("no integer identifier found in '%s'\n",readbuffer.c_str());
    readbuffer.erase(0, end_pos+1);

    // extract flag
    start_pos = readbuffer.find('\t', end_pos+1) + 1;
    end_pos = readbuffer.find('\t', start_pos);
    if(end_pos == string::npos)
	Rf_error("failed to find sam flag in '%s'\n", readbuffer.c_str());
    readflag = atoi(readbuffer.substr(start_pos, end_pos-start_pos).c_str());

    // set readIsMapped and readIsPaired
    readIsMapped = !(readflag & BAM_FUNMAP);
    readIsPaired = (readflag & BAM_FPAIRED);
    readIsMapped2 = !(readflag & BAM_FMUNMAP);

    // read second if paired and mate is mapped
    if((readIsMapped && readIsPaired && readIsMapped2) ||
       (!readIsMapped && readIsPaired && !readIsMapped2)) {
	// read line
	getline (fh, readbuffer2, '\n');
	if(fh.eof() || !fh.good())
	    Rf_error("error reading second alignment of pair from %s\n", fname);

        // remove \r if exists (for windows)
        if(readbuffer2[readbuffer2.size()-1] == '\r')
            readbuffer2.erase(readbuffer2.size()-1, 1);

	// extract id
	end_pos = readbuffer2.find('_');
	if(end_pos != string::npos)
	    readid2 = atoi(readbuffer2.substr(0, end_pos).c_str());
	else
	    Rf_error("no integer identifier found in '%s'\n",readbuffer2.c_str());
	readbuffer2.erase(0, end_pos+1);

	// extract flag
	start_pos = readbuffer2.find('\t', end_pos+1) + 1;
	end_pos = readbuffer2.find('\t', start_pos);
	if(end_pos == string::npos)
	    Rf_error("failed to find sam flag in '%s'\n", readbuffer2.c_str());
	readflag2 = atoi(readbuffer2.substr(start_pos, end_pos-start_pos).c_str());

	// check if paired
	if(readid!=readid2 || !(readflag2 & BAM_FPAIRED)) {
	    Rf_error("unexpected alignment when reading second of a pair\n");

	} else {
	    // adjust readIsMapped
	    readIsMapped = (readIsMapped || readIsMapped2);
	}

    } else {
	readbuffer2.clear();
    }

    return 0;
}

// read and store alignments from fh until readid == id, and then until readid != id
int samFile::advance(int id) {
    //cout << "advancing(" << id << "), top: " << (queue.empty() ? -1 : queue.top().id) << endl;

    static int nr;
    static streampos filePos;

    if(!fh.eof() && (queue.empty() || queue.top().id != id)) {
	// do nothing if EOF reached or id is already on queue.top()

	readid = -1;
	nr = 0;
	while (readid != id) {
	    // read next alignment
	    if(this->getNextAln())
		break;

	    // store in queue
	    queue.push(idLine(readid, readIsMapped, readbuffer, readbuffer2));
	    //cout << "\tjust stored " << readid << endl;
	    nr++;
	}

	// read all alignments with that id
	while (readid == id) {
	    // read next alignment
	    filePos = fh.tellg();
	    if(this->getNextAln())
		break;

	    if(readid == id) { // same id
		// store in queue
		queue.push(idLine(readid, readIsMapped, readbuffer, readbuffer2));
		//cout << "\tjust stored " << readid << endl;
		nr++;

	    } else {
		// next id found; undo last getline
		//cout << "\tjust ignored " << readid << endl;
		fh.seekg(filePos);
		if(fh.fail() || fh.bad())
		    Rf_error("failed to seek to new position in sam file");
	    }
	}
    } else {
	nr = 0;
    }

    //cout << "\t" << nr << " alignments parsed, new top: " << (queue.empty() ? -1 : queue.top().id) << endl;
    return (int)(queue.size());
}

// output all stored alignments for 'id' and remove them from memory, store unmapped reads in 'unmapped'
int samFile::flush_simple(int id, ofstream &outfh, map<int, string> &unmapped) {
    static int numberFlushed;

    numberFlushed = 0;
    while(!queue.empty() && queue.top().id == id) {

	// output
	if(queue.top().isMapped) {
	    // ... to file
	    outfh << queue.top().line << '\n';
	    if(! queue.top().line2.empty())
		outfh << queue.top().line2 << '\n';
	    numberFlushed++;

	} else if(unmapped.count(id) == 0) {
	    // ... to 'unmapped' map for later nonredundant output
	    if(! queue.top().line2.empty())
		unmapped[id] = queue.top().line + '\n' + queue.top().line2;
	    else
		unmapped[id] = queue.top().line;
	}

	// remove from queue
	queue.pop();
    }

    return numberFlushed;
}

/* bisulfite-version of samFile::flush_simple, this does in addition:
   - correct identifier (remove readseq_ prefix)
   - fix sequence (replace read sequence with the one contained in the id)
   - fix MD tag (remove MD tag if present, as it could be inaccurate due to the read sequence exchange)

   it is assumed that the library was --fr and has been changed to --ff for alignment, therefore:
   - for paired alignments, first reads: change strand of next fragment in flag
   - for paired alignments, second reads: change strand of fragment in flag, reverse-complement the sequence  */
int flush_bisulfite(int id, ofstream &outfh, map<int, string> &unmapped, vector<idLine> &mapped, bool addId) { 
    static int numberFlushed;
    numberFlushed = 0;
    static idLine currenttop;
    static vector<idLine>::size_type i, count;
    i = 0;
    count = mapped.size();

    while(i < count){
	if(addId){
	    currenttop = mapped[i];
	    i++;
	} else {
	    currenttop = mapped[(unsigned long)(unif_rand()*count)];
	    i = count;
	}

	_fix_FLAGs_and_sequences(currenttop);
	
	// output
	if(currenttop.isMapped) {
	    // ... to file
	    if(addId) {
		outfh << id << '_' << currenttop.line << '\n';
		if(! currenttop.line2.empty())
		    outfh << id << '_' << currenttop.line2 << '\n';
	    } else {
		outfh << currenttop.line << '\n';
		if(! currenttop.line2.empty())
		    outfh << currenttop.line2 << '\n';
	    }
	    numberFlushed++;	
	}
    }

    return numberFlushed;
}

// allele-specific-version of 'flush_simple', this does in addition:
// - add 'allele-tag'
int flush_allele(int id, ofstream &outfh, map<int, string> &unmapped, idLine &currenttop, char tag) {
    static int numberFlushed;
    numberFlushed = 0;

    // add allele-tag output to file
    if(!currenttop.line2.empty()){
	outfh << currenttop.line << '\t' << "XV:A:" << tag << '\n';
	outfh << currenttop.line2 << '\t' << "XV:A:" << tag << '\n';
    } else
	outfh << currenttop.line << '\t' << "XV:A:" << tag << '\n';
    numberFlushed++;

    return numberFlushed;
}

// output unmapped for 'id' and remove from memory (just once per read)
// NOTE: The parameter 'id' is unused currently. The parameter 'unmapped' contains 
//       always only one element, which correspond to the current id.
//       The 'id' parameter exist because of historical reason.
//       It could be removed but we keep it in case we change the algorithm.
int samFile::flush_unmapped(int id, ofstream &outfh, map<int, string> &unmapped, int n) {
    static map<int, string>::iterator it;
    static int numberFlushed;
    numberFlushed = 0;

    if(n==0) {
	numberFlushed = (int)(unmapped.size());

	for(it=unmapped.begin(); it != unmapped.end(); it++)
	    outfh << it->second << '\n';
    }
    unmapped.clear();

    return numberFlushed;
}

int samFile::get_alignments_bisulfite(int id, int bisQueue, vector<idLine> &mapped, map<int, string> &unmapped, bool output, bool addId) {
    static idLine currenttop;
    static char idbuffer[64];

    while(!queue.empty() && queue.top().id == id) {
	currenttop = queue.top();

	// check if unmapped
	if(!currenttop.isMapped){
	    if(unmapped.count(id) == 0) {
		// add to 'unmapped' map for later nonredundant output
		_replace_sequence(currenttop.line, false);
		if(!currenttop.line2.empty())
		    _replace_sequence(currenttop.line2, false);
		if(addId) {
		    sprintf(idbuffer, "%i", id);
		    if(! currenttop.line2.empty())
			unmapped[id] = ((string)idbuffer + '_' + currenttop.line + '\n' + idbuffer + '_' + currenttop.line2);
		    else
			unmapped[id] = ((string)idbuffer + '_' + currenttop.line);
		} else {
		    if(! currenttop.line2.empty())
			unmapped[id] = (currenttop.line + '\n' + currenttop.line2);
		    else
			unmapped[id] = (currenttop.line);
		}
	    }
	} else if(output){ // if output == true then ...
	    // ... set the bisulfite queue member variable
	    currenttop.bisQueue = bisQueue;
	    // ... and add mapped alignments to 'mapped' vector
	    mapped.push_back(currenttop);
	}

	// remove from queue
	queue.pop();
    }
    return EXIT_SUCCESS;
}

int samFile::get_alignments_allele(int id, vector<idLine> &mapped, map<int, string> &unmapped) {
    static idLine currenttop;

    while(!queue.empty() && queue.top().id == id) {
	currenttop = queue.top();

	// check if unmapped
	if(!currenttop.isMapped){
	    if(unmapped.count(id) == 0) {
		// add to 'unmapped' map for later nonredundant output
		if(! currenttop.line2.empty())
		    unmapped[id] = (currenttop.line + '\n' + currenttop.line2);
		else
		    unmapped[id] = (currenttop.line);
	    }
	} else {
	    // ... add mapped alignments to 'mapped' vector
	    mapped.push_back(currenttop);
	}

	// remove from queue
	queue.pop();
    }
    return EXIT_SUCCESS;
}

int samFile::get_nm_tag(int id){
    if(queue.top().id == id && queue.top().isMapped) {
	// mapped get edit distance
	return _get_nm_tag(queue.top());
    } else {
	return MAX_NM;
    }
}

int _get_nm_tag(const idLine &alignment){
    static int nm;
    nm = MAX_NM;
    static size_t pos;

    // get edit distance
    pos = alignment.line.find("NM:i:")+5;
    nm = alignment.line[pos] - '0';
    // if paired then get edit distance of pair and add up
    if(!alignment.line2.empty()){
	pos = alignment.line2.find("NM:i:")+5;
	nm = nm + (alignment.line2[pos] - '0');
    }

    return nm;
}

inline char complement(char element) {
    static const char charMap[] = {
	'T', 'V', 'G', 'H', 'N', 'N', 'C', 'D', 'N', 'N', 'M', 'N', 'K',
	'N', 'N', 'N', 'N', 'Y', 'S', 'A', 'A', 'B', 'W', 'N', 'R', 'N'
    };
    return charMap[toupper(element) - 'A'];
}

void _reverse_complement(string &seq) {
    reverse(seq.begin(), seq.end());
    transform(seq.begin(), seq.end(), seq.begin(), complement);
}

void _replace_sequence(string &line, bool revcomp) {
    static size_t start_pos, end_pos;
    static int i;
    static string origseq;

    // get sequence from beginning of line (anyting up to first '_') and store in origseq
    end_pos = line.find('_');
    if(end_pos != string::npos)
	origseq = line.substr(0, end_pos);
    else
	Rf_error("no read sequence found in '%s'\n",line.c_str());
    line.erase(0, end_pos+1);

    // reverse-complement?
    if(revcomp)
	_reverse_complement(origseq);

    // replace seq (10th tab-delimited field in line) by the one in origseq
    start_pos = line.find('\t') + 1;
    for(i=0; i<8; i++)
	start_pos = line.find('\t', start_pos) + 1;
    end_pos = line.find('\t', start_pos);
    if(start_pos != string::npos && end_pos != string::npos)
	line.replace(start_pos, end_pos-start_pos, origseq);
    else
	Rf_error("error finding sequence column in '%s'\n",line.c_str());
}

void _remove_MD_tag(string &line) {
    static size_t start_pos, end_pos;

    start_pos = line.rfind("MD:Z:", string::npos, 5);
    if(start_pos != string::npos && line[start_pos-1] == '\t') {
	end_pos = line.find('\t', start_pos+5);
	if(end_pos != string::npos)
	    line.erase(start_pos, end_pos-start_pos+1);
	else
	    line.erase(start_pos, line.length()-start_pos);
    }
}

void _fix_FLAGs_and_sequences(idLine &currenttop) {
    static bool revcomp;
    static char tagbuffer[64];

    revcomp = currenttop.bisQueue % 2 ? true : false;

    _replace_sequence(currenttop.line, revcomp);

    sprintf(tagbuffer, "\tXQ:i:%i", currenttop.bisQueue);
    currenttop.line += tagbuffer;

    if(! currenttop.line2.empty()) { 
	currenttop.line2 += tagbuffer;
	_replace_sequence(currenttop.line2, revcomp);
    }
}

// QUESTION return size_t?
int _fix_identical_locus(vector<idLine> &mapped){
    static map<string,int> locus;
    locus.clear();
    static map<string,int>::iterator locus_it;
    static size_t start_pos, end_pos;
    static vector<idLine>::size_type i;
    static int j, comp;
    static string rname, pos, rname2, pos2;
    string key;
    static bool rm_first = false; // remove first of identical locus 
    idLine curr;

    // get locus to create the key
    for(i=0; i < mapped.size(); i++){
	// get reference name
	start_pos = 0;
	for(j=0; j<2; j++)
	    start_pos = mapped[i].line.find('\t', start_pos) + 1;
	end_pos = mapped[i].line.find('\t', start_pos);
	rname = mapped[i].line.substr(start_pos, end_pos-start_pos);
	// get position
	start_pos = end_pos + 1;
	end_pos = mapped[i].line.find('\t', start_pos);
	pos = mapped[i].line.substr(start_pos, end_pos-start_pos);
	// check if paired
	if(mapped[i].line2.empty()){
	    // add key to map
	   key = rname+pos;
	} else {
	    // get reference name of mate
	    start_pos = 0;
	    for(j=0; j<2; j++)
		start_pos = mapped[i].line2.find('\t', start_pos) + 1;
	    end_pos = mapped[i].line2.find('\t', start_pos);
	    rname2 = mapped[i].line2.substr(start_pos, end_pos-start_pos);
	    // get position of mate
	    start_pos = end_pos + 1;
	    end_pos = mapped[i].line2.find('\t', start_pos);
	    pos2 = mapped[i].line2.substr(start_pos, end_pos-start_pos);
	    // sort location and add key to map
	    comp = rname.compare(rname2);
	    if(comp == 0){
		if(pos.compare(pos2) <= 0)
		    key = rname+pos+rname2+pos2;
		else 
		    key = rname2+pos2+rname+pos;
	    } else if(comp < 0){
		key = rname+pos+rname2+pos2;
	    } else {
		key = rname2+pos2+rname+pos;
	    }
	}
	locus_it = locus.find(key);
      	if(locus_it == locus.end()){
	    locus.insert(pair<string,int>(key,(const int)i));
       	} else {
	    // remove one of the duplicated location
	    if(rm_first){
		curr = mapped[i];
		mapped.insert(mapped.begin()+locus_it->second, curr);
		mapped.erase(mapped.begin()+locus_it->second+1);
	    }
	    mapped.erase(mapped.begin()+(const long)i);
	    rm_first = !rm_first;
	    i--;
	}
    }

    // return number of identical locus
    return (int)locus.size();
}

int _make_unmapped_alignment(int id, idLine &mapped_alignment, map<int, string> &unmapped, bool addId, bool replaceSeq) {
    string qname, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual;
    static char int_buffer[64];
    int flag;
    static size_t end_pos;
    static bool reverse;
    string line, line2;
    static string str_buffer;

    // create unmapped read
    istringstream iss(mapped_alignment.line);
    iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
    // modify flag
    reverse = flag & BAM_FREVERSE;
    if(mapped_alignment.line2.empty())
	flag = BAM_FUNMAP + (flag & ~(BAM_FPROPER_PAIR + BAM_FUNMAP + BAM_FMUNMAP + BAM_FREVERSE + BAM_FMREVERSE));
    else
	flag = BAM_FUNMAP + BAM_FMUNMAP + (flag & ~(BAM_FPROPER_PAIR + BAM_FUNMAP + BAM_FMUNMAP + BAM_FREVERSE + BAM_FMREVERSE));
    sprintf(int_buffer, "%i", flag);
    // replace sequence
    if(replaceSeq){
	// bisulfite mode: no reverse necessary (taking seq from qname)
	// get sequence from beginning of line (anyting up to first '_') and store in origseq
	end_pos = qname.find('_');
	if(end_pos != string::npos){
	    seq = qname.substr(0, end_pos);
	    qname.erase(0, end_pos+1);
	} else
	    Rf_error("no read sequence found in '%s'\n", qname.c_str());
    } else if(reverse){
	// allelic mode
	_reverse_complement(seq);
    }
    // reverse qualitiy
    if(reverse)
	std::reverse(qual.begin(), qual.end());
    line = qname + "\t" + int_buffer + "\t*\t0\t0\t*\t*\t0\t0\t" + seq + "\t" + qual;

    // create unmapped mate read
    if(! mapped_alignment.line2.empty()){
	istringstream iss(mapped_alignment.line2);
	iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
	// modify flag
	reverse = flag & BAM_FREVERSE;
	flag = BAM_FUNMAP + BAM_FMUNMAP + (flag & ~(BAM_FPROPER_PAIR + BAM_FUNMAP + BAM_FMUNMAP + BAM_FREVERSE + BAM_FMREVERSE));
	sprintf(int_buffer, "%i", flag);
	// replace sequence
	if(replaceSeq){
	    // bisulfite mode: no reverse necessary (taking seq from qname)
	    // get sequence from beginning of line (anyting up to first '_') and store in origseq
	    end_pos = qname.find('_');
	    if(end_pos != string::npos){
		seq = qname.substr(0, end_pos);
		qname.erase(0, end_pos+1);
	    } else
		Rf_error("no read sequence found in '%s'\n", qname.c_str());
	} else if(reverse){
	    // allelic mode
	    _reverse_complement(seq);
	}
	// reverse qualitiy
	if(reverse)
	    std::reverse(qual.begin(), qual.end());
	line2 = qname + "\t" +  int_buffer + "\t*\t0\t0\t*\t*\t0\t0\t" + seq + "\t" + qual;

	// check order of read1 and read2. Read1 should be first in the output.
	if(flag & BAM_FREAD1){
	    str_buffer = line;
	    line = line2;
	    line2 = str_buffer;
	}
    }

    if(addId) {
	sprintf(int_buffer, "%i", id);
	if(! mapped_alignment.line2.empty())
	    unmapped[id] = ((string)int_buffer + '_' + line + '\n' + int_buffer + '_' + line2);
	else
	    unmapped[id] = ((string)int_buffer + '_' + line);
    } else {
	if(! mapped_alignment.line2.empty())
	    unmapped[id] = (line + '\n' + line2);
	else
	    unmapped[id] = (line);
    }

    return EXIT_SUCCESS;
}

int _fix_half_mapper(vector<idLine> &mapped, map<int, string> &unmapped){
    static size_t start_pos, end_pos;
    static vector<idLine>::size_type i;
    static int id, flag = 0;
    string line1, line2;
    vector<idLine>::iterator it = mapped.begin();

     for(i=0; i < mapped.size(); i++){
	// get flag
	start_pos = mapped[i].line.find('\t');
	end_pos = mapped[i].line.find('\t', start_pos+1);
	if(start_pos == string::npos || end_pos == string::npos)
	    Rf_error("failed to find sam flag in '%s'\n", mapped[i].line.c_str());
	flag = atoi(mapped[i].line.substr(start_pos, end_pos-start_pos).c_str());
	if((flag & BAM_FPAIRED) && (flag & BAM_FMUNMAP) && mapped[i].line2.empty()){
	    // half mapper
	    id = mapped[i].id;
	    if(flag & BAM_FREAD2)
		line2 = mapped[i].line;
	    else
		line1 = mapped[i].line;
	    // remove
	    mapped.erase(it+(const long)i);
	    i--;
	}
     }
     // make unmapped read from half mapper
     if(!line1.empty()){
	 idLine halfmapper = idLine(id, true, line1, line2);
	 _make_unmapped_alignment(id, halfmapper, unmapped, false, false);
     }

    return EXIT_SUCCESS;
}

// copy SAM header from samfile and write to output file handle
int _copy_header(const char *headerin, ofstream &outfh) {
    ifstream fh;
    string linebuffer;

    fh.open(headerin, ifstream::in);
    if(! fh.good()) {
	Rf_error("error opening file '%s'\n",headerin);
    } else {
	// copy header
	while(fh.good()) {
	    getline (fh, linebuffer, '\n');
	    if( linebuffer[0] == '@' )
		outfh << linebuffer << '\n';
	    else
		break;
	}
    }
    fh.close();
    return 0;
}

// define functions for output generation (will be called through pointers according to the merge mode, defined by 'bisulfiteMode' and 'alleleMode'
int writeOutput_simple(int id, samFile **samf, int nsamf, ofstream &outfh, map<int, string> &unmapped, int maxhits) {
    static int n, i;
    n = 0;
    for(i=0; i<nsamf; i++)
	n += samf[i]->flush_simple(id, outfh, unmapped);
    return n;
}

int writeOutput_bisulfite_core(int id, samFile **samf, int nsamf, ofstream &outfh, map<int, string> &unmapped, int maxhits, bool addId) {
    static int n, i, min_nm, curr_nm, count;
    vector<idLine> mapped;
    n = 0; // number of flashed alignments
    curr_nm = MAX_NM; // current nm tag value
    min_nm = MAX_NM; // smallest nm tag value

    // read nm tag and find smallest value
    for(i=0; i<nsamf; i++){
	curr_nm = samf[i]->get_nm_tag(id);
	if(curr_nm < min_nm){
	    // better mapped alignments found empty vector 
	    min_nm = curr_nm;
	    mapped.clear();
	    samf[i]->get_alignments_bisulfite(id, i, mapped, unmapped, true, addId);	    
	} else if(curr_nm == min_nm){
	    // equal mapped alignments add them to the vector 
	    samf[i]->get_alignments_bisulfite(id, i, mapped, unmapped, true, addId);
	} else {
	    // worse mapped alignments to not output
	    samf[i]->get_alignments_bisulfite(id, i, mapped, unmapped, false, addId);
	}
    }

    // if undirected bisulfite fix alignments with identical locus
    if(nsamf > 2)
	_fix_identical_locus(mapped);

    // fix halfmapper not needed for bisulfit (no halfmapper generated by bowtie1)

    count =  (int)mapped.size();    
    // if there are mapped alignments then output 
    if(count > 0){
	if(!addId && count > maxhits){ // if addId = true then before allele. output all best mapped read
	    if(unmapped.count(id) == 0){
		// more than maxhits mapped alignments -> output as unmapped
		_make_unmapped_alignment(id, mapped[0], unmapped, addId, true);
	    }
	} else {
	    // equal or less than maxhits alignments -> output
	    n += flush_bisulfite(id, outfh, unmapped, mapped, addId);
	}
    }

    return n;
}

int writeOutput_bisulfite(int id, samFile **samf, int nsamf, ofstream &outfh, map<int, string> &unmapped, int maxhits) {
    static int n = 0;
    n = writeOutput_bisulfite_core(id, samf, nsamf, outfh, unmapped, maxhits, false);
    return n;
}

int writeOutput_bisulfite_before_allele(int id, samFile **samf, int nsamf, ofstream &outfh, map<int, string> &unmapped, int maxhits) {
    static int n = 0;
    n = writeOutput_bisulfite_core(id, samf, nsamf, outfh, unmapped, maxhits, true);
    return n;
}

int writeOutput_allele(int id, samFile **samf, int nsamf, ofstream &outfh, map<int, string> &unmapped, int maxhits) {
    if(nsamf != 2)
	Rf_error("Only two input files are allowed for allele specific mode.");

    static int n; // number of alignment writen to the output
    n = 0;
    static bool allele = true; // alternate allele, boolean which input file should be choosen
    static int nmR, nmA, countR, countA;

    // get alignments from the queue
    vector<idLine> mappedR;
    vector<idLine> mappedA;
    samf[0]->get_alignments_allele(id, mappedR, unmapped);
    samf[1]->get_alignments_allele(id, mappedA, unmapped);

    // fix half mapper
    _fix_half_mapper(mappedR, unmapped);
    _fix_half_mapper(mappedA, unmapped);

    // find alignment with fewest mismatch and add allele tag
    countR =  (int)mappedR.size();
    countA =  (int)mappedA.size();
    nmR = MAX_NM;
    nmA = MAX_NM;
    if(countR > 0)
	nmR = _get_nm_tag(mappedR[0]);
    if(countA > 0)
	nmA = _get_nm_tag(mappedA[0]);

    if(nmR != nmA){
	if(nmR < nmA){
	    // ref less mismatch
	    if(countR > maxhits) // over mapped
		_make_unmapped_alignment(id, mappedR[0], unmapped, false, false);
	    else
		n += flush_allele(id, outfh, unmapped, mappedR[(unsigned long)(unif_rand()*countR)], 'R');
	} else{
	    // alternate less mismatch
	    if(countA > maxhits) // over mapped
		_make_unmapped_alignment(id, mappedA[0], unmapped, false, false);
	    else
		n += flush_allele(id, outfh, unmapped, mappedA[(unsigned long)(unif_rand()*countA)], 'A');
	}
    } else {
	// both same number of mismatch
	if(countR > maxhits || countA > maxhits) // over mapped
	    _make_unmapped_alignment(id, mappedR[0], unmapped, false, false);
	else if(countR > 0 && countA > 0){
	    if(allele)
		n += flush_allele(id, outfh, unmapped, mappedR[(unsigned long)(unif_rand()*countR)], 'U');
	    else
		n += flush_allele(id, outfh, unmapped, mappedA[(unsigned long)(unif_rand()*countA)], 'U');
	    allele = !allele; // switch allele 
	}
    }
 
    return n;
}

int _merge_reorder_sam(const char** fnin, int nin, const char* fnout, int mode, int maxhits) {
    int i, id, n, e, currQueueSize = 0, maxQueueSize = 0;
    map <int, string> unmapped;

    typedef int (*OUTPUTFUNCTION) (int, samFile**, int, ofstream&, map<int, string>&, int); // pointer to the output function
    OUTPUTFUNCTION writeOutput = NULL;

    // set function pointer for output according to 'mode'
    switch (mode) {
    case 0: // 0 : simple (any number of input files, remove id from QNAME)
	writeOutput = writeOutput_simple; break;
    case 1: // 1 : bisulfite-mode (2 or 4 input files, remove id and seq from QNAME, replace SEQ, remove MD-tag)
	writeOutput = writeOutput_bisulfite; break;
    case 2: // 2 : allele-specific-mode (2 input files, remove id, add 'allele'-tag)
	writeOutput = writeOutput_allele; break;
    case 3: // 3 : bisulfite-mode, to be followed by allele-specific-mode (as '1', but leave id in QNAME)
	writeOutput = writeOutput_bisulfite_before_allele; break;
    default:
	Rf_error("'mode' must be 0, 1, 2 or 3");
    }

    // open output file
    ofstream outfile (fnout, ofstream::out | ofstream::binary);
    if(! outfile.good())
	Rf_error("error opening output file: %s\n", fnout);

    // copy header from first input file
    if(_copy_header(fnin[0], outfile))
	Rf_error("error copying header from %s\n", fnin[0]);
    
    // open sam files
    samFile **samfiles = new samFile*[nin];
    for(i=0; i<nin; i++)
	samfiles[i] = new samFile(fnin[i]);

    // main loop over identifiers (1...allEof)
    id = 1;
    while (!samFile::allEof()) {
	// forward input files until current identifier is found
	// unmapped reads from all samfiles are collected in unmappedQueue
	for(i=0; i<nin; i++)
	    currQueueSize = samfiles[i]->advance(id);
	//cout << samfiles[i]->advance(id) << endl;
	if(currQueueSize > maxQueueSize)
	    maxQueueSize = currQueueSize;

	// output alignments for current identifier
	//   and remove alignments for current identifier from memory
	n = writeOutput(id, samfiles, nin, outfile, unmapped, maxhits);

	// output or delete unmapped
	samFile::flush_unmapped(id, outfile, unmapped, n);
	
	// increase current identifier
	id++;
    }

    // files are done; flush memory
    e = 0;
    for(i=0; i<nin; i++)
	if(!samfiles[i]->isEmpty())
	    e++;
    //cout << "STILL HAVE " << e << endl;
    while(e>0) {
	// output alignments for current identifier
	//   and remove alignments for current identifier from memory
	n = writeOutput(id, samfiles, nin, outfile, unmapped, maxhits);
	for(i=0; i<nin; i++) {
	    if(samfiles[i]->isEmpty())
		e--;
	}

	// output or delete unmapped
	samFile::flush_unmapped(id, outfile, unmapped, n);
	
	// increase current identifier
	id++;
    }

    // close sam files and clean up
    for(i=0; i<nin; i++)
        delete samfiles[i];
    delete[] samfiles;

    return maxQueueSize;
}


#ifdef __cplusplus
extern "C" {
#endif

SEXP merge_reorder_sam(SEXP infiles, SEXP outfile, SEXP mode, SEXP maxhits) {
    if (!Rf_isString(infiles))
	Rf_error("'infiles' must be a character vector");
    if (!Rf_isString(outfile) || 1 != Rf_length(outfile))
	Rf_error("'outfile' must be a single character value");
    if (!Rf_isInteger(mode) || 1 != Rf_length(mode))
        Rf_error("'mode' must be integer(1)");
    if (!Rf_isInteger(maxhits) || 1 != Rf_length(maxhits))
        Rf_error("'maxhits' must be integer(1)");

    int i = 0, res = 0, nbIn = Rf_length(infiles), mode_int = Rf_asInteger(mode);

    if (mode_int < 0 || mode_int > 3)
        Rf_error("'mode' must be 0, 1, 2 or 3");
    if ((mode_int == 1 || mode_int == 3) && (nbIn != 2 && nbIn != 4))
        Rf_error("in mode=1 and mode=3 (bisulfite), there must be exactly 2 or 4 input files");
    if (mode_int == 2 && nbIn != 2)
        Rf_error("in mode=2 (allele-specific), there must be exactly 2 input files");

    const char **inf = (const char**) R_Calloc(Rf_length(infiles), char*);
    for(i=0; i<nbIn; i++)
	inf[i] = Rf_translateChar(STRING_ELT(infiles, i));
    GetRNGstate(); // prepare use of R random number generator
    res = _merge_reorder_sam(inf, nbIn, Rf_translateChar(STRING_ELT(outfile, 0)), mode_int, Rf_asInteger(maxhits));
    PutRNGstate(); // finish use of R random number generator
    R_Free(inf);

    return Rf_ScalarInteger(res);
}

#ifdef __cplusplus
}
#endif

/*
int main(int argc, const char **argv) {
    _merge_reorder_sam(argv+1, argc-3, argv[argc-2], atoi(argv[argc-1]));
}
*/

