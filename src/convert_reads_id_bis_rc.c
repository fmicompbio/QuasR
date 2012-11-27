#include "convert_reads_id_bis_rc.h"
#include "utilities.h"

KSEQ_INIT(gzFile, gzread)

/*! @function
  @abstract  Convert fasta or fastq sequence to a three letter base space and add the original sequence to the id.
  @param  inFile  character(1) Input name of the fasta/fastq file
  @param  outFile  character(1) Output name of the fasta/fastq file
  @param  fromToCharacter  character(2) or NULL; two single capital letter character of the base from which to which should be converted. If NULL no conversion is done.
  @param  reverseComplement  logical(1) if TRUE the sequence is transformed to the reversecomplement before the conversion 
 */

SEXP convert_reads_id_bis_rc(SEXP inFile, SEXP outFile, SEXP fromToCharacter, SEXP reverseComplement)
{
    int convert = 1;
    // check parameters
    if(!Rf_isString(inFile) || 1 != Rf_length(inFile)){
        Rf_error("'infile' must be character(1)");
    }
    if(!Rf_isString(outFile) || 1 != Rf_length(outFile)){
        Rf_error("'outfile' must be character(1)");
    }
    if(Rf_isNull(fromToCharacter)){
        convert = 0;
    } else if(!Rf_isString(fromToCharacter) || 2 != Rf_length(fromToCharacter)){
        Rf_error("'fromToCharacter' must be character(2)");
    }
    if(!Rf_isLogical(reverseComplement) || 1 != Rf_length(reverseComplement)){
        Rf_error("'reverseComplement' must be logical(1)");
    }
    
    char from = 0, to = 0; 
    if(convert){
	from = Rf_translateChar(STRING_ELT(fromToCharacter, 0))[0];
        to = Rf_translateChar(STRING_ELT(fromToCharacter, 1))[0];
	if(64 > from || from > 91)
	    Rf_error("'from' must be a capital letter");
	if(64 > to || to > 91)
	    Rf_error("'to' must be a capital letter");
    }

    gzFile fin;
    fin = gzopen(Rf_translateChar(STRING_ELT(inFile, 0)), "r");
    int l;
    kseq_t *seq = kseq_init(fin);
        
    FILE * fout = NULL;
    fout = fopen(Rf_translateChar(STRING_ELT(outFile, 0)), "w");
    if(fout == NULL)
        Rf_error("Error in opening the output file  %s", Rf_translateChar(STRING_ELT(outFile, 0)));

    uint64_t count = 1;
    char first_char = '>';
    int offset = (int)'A';
    char translate[] = "ANCNNNGNNNNNNNNNNNNTNNNNNNNNNNNNaNcNNNgNNNNNNNNNNNNtNNNNNN";
    if(convert){
	// create translation array
	translate[(int)from - offset] = to; // convert capital letter
	translate[(int)from - offset + 'a' -offset] = to + 'a' -offset; // convert small letter
    }

    while ((l = kseq_read(seq)) >= 0){
        // define first character of the seq name line
        if(seq->qual.l)
            first_char = '@';
        else
            first_char = '>';

        // remove \r if exists (for windows)
        if(seq->seq.s[seq->seq.l-1] == '\r'){
            seq->seq.s[seq->seq.l-1] = '\0';
            seq->seq.l = seq->seq.l-1;
            if(seq->qual.l){
                seq->qual.s[seq->qual.l-1] = '\0';
                seq->qual.l = seq->qual.l-1;
            }
        }

        // reverse complement seq and reverse quality
        if(Rf_asLogical(reverseComplement)){
            _complement(seq->seq.s, seq->seq.l);
            _reverse(seq->seq.s, seq->seq.l);
            _reverse(seq->qual.s, seq->qual.l);
        }
        
	// write seq-name and convert sequence
	if(convert){
	    // write new seq-name in format 'count_unconvertedSequence_seqname'
	    if (seq->comment.l)
		fprintf(fout, "%c%"PRIu64"_%s_%s %s\n", first_char, count, seq->seq.s, seq->name.s, seq->comment.s);
	    else
		fprintf(fout, "%c%"PRIu64"_%s_%s\n", first_char, count, seq->seq.s, seq->name.s);
	    // convert sequence
	    for(int i = 0; i < seq->seq.l; i++){
		seq->seq.s[i] = translate[(int)seq->seq.s[i] - offset];
	    }	    
	} else {
	    // write new seq name in format 'count_seqname'
	    if (seq->comment.l)
		fprintf(fout, "%c%"PRIu64"_%s %s\n", first_char, count, seq->name.s, seq->comment.s);
	    else
		fprintf(fout, "%c%"PRIu64"_%s\n", first_char, count, seq->name.s);
	}
        // write sequence
        fprintf(fout, "%s\n", seq->seq.s);
        // write qualtity
        if (seq->qual.l)
            fprintf(fout, "+\n%s\n", seq->qual.s);
        count++;
     }

    kseq_destroy(seq);
    gzclose(fin);
    fclose(fout);
    return Rf_ScalarInteger(0);
} 
