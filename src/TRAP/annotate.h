// ============================================================================
// This code was taken and slightly adapted from TRAP/ANNOTATE v3.04
// http://trap.molgen.mpg.de/cgi-bin/home.cgi
// Not original code of PROmiRNA
// ============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include "annotate_sequence.h"
#include "matrix.h"
#include "score_calc.h"
#include <time.h>
#include <cmath>
#include <cstdlib>
#include "string.h"
#include <algorithm>

namespace trap // NAMESPACE BECAUSE OF GLOBAL VARIABLES
{

using namespace std;

//function declarations for this class
void initialize_Sequenceinput_GFFoutput();
void validate_options();
void annotate(vector<double>, const char *, const char *, float);

//variables
char gap_char;
double gc_content; //external gc_content variable
double score_thresh=0; //score cutoff hardwired default value '0'
double pscount=1; //pseudo count hard wired to 1
string thresh_type="balanced"; //for the 3-types cutoff, default is balanced
double fact_lambda=0.7; //lambda in affinity based annotation hardwired to '0.7'
double fact_lnR0;
long region_size;
int number_of_top_scores;
string mask_handle = "average";
string single_strand_direct;
vector<double> GC_vector;
vector<class motifMatrix> pssm_motif_matrices;
vector<class motifMatrix> psem_motif_matrices;

//vector that contains positions of the motif to discard
vector<int> discard_pos;
string discard_pos_param;

// flags serving different purposes
bool GAP=false;
bool DEF_GC=false;
bool GC_RANGE = false;
bool NAME_FLAG=false;
bool ANNOTATION_AVERAGE=false;
bool MASK_SCORING=true;
bool DISCARD_POS=false;
bool SINGLE_STRAND = false;
bool SS_FORWARD = true;
bool AFFINITY_SCORE=false;
bool EXTN_GC_FLAG=false;
bool LARGEST=false;
bool PSSM=false;
bool PSCM=false;
bool PSEM=false;
bool EXIT=false;
bool USER_SHIFT=false;
bool EXTN_SEQ_FILE=false;
bool EXTN_GFF_FILE=false;
bool FACT_LAMBDA=false;
bool FACT_LNR0=false;
bool EXTN_THRESH=false;
bool THRESH_TYPE = false;
bool DEBUG=false;
bool REGION=false;
bool TABLE=false;
int region_shift=0;

char *matFile,*seqFile,*ps_count,*regSize,*regShift,*nTopScore; //*op_file
vector<string> name_flags;
istream *user_seq_file; //stream object for input of sequence
ostream *user_gff_file;	//stream object for output of gff

vector<double> tmp_affinities;

///////////////-----------------------------------------------------/////////////////////////
// Function for I/O (adapted)
void initialize_Sequenceinput_GFFoutput()
{

  user_seq_file=&cin;  // initializing the stream object for reading sequences with STANDARD INPUT;
  if(EXTN_SEQ_FILE){  // if external paramater -s is specified... override the stream object with

    user_seq_file=new ifstream(seqFile);   // the file for the input file specified by the user

    // check if the file exists
    if(user_seq_file->fail()){
      cerr<<"ERROR: SEQUENCE FILE NOT FOUND"<<endl;
      exit(1);
    }

  }

  user_gff_file=&cout;	    // initializing the stream object for gff output as STANDARD OUTPUT
  // if(EXTN_GFF_FILE)	    // if external paramater is specified... override the stream object with
  // user_gff_file=new ofstream(op_file);	     // the file for the output file specified by the user

  // start the output file with the command used
  // *user_gff_file <<"#command_line= ";
  // for (int la=0; la<argc; la++) *user_gff_file << argv[la] << " "; *user_gff_file<<endl;

  // headers for the columns of the output
  // if (!TABLE){
  // 	  *user_gff_file <<"#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes"; ; *user_gff_file<<endl;
  // }
}

///////////////-----------------------------------------------------/////////////////////////
// Routine to check for the validity of the options of the tool
// NOTE: Commented code is only needed if more options are added for PROmiRNA (anyway the code will be kept)

void validate_options()
{

// 	if(!EXTN_SEQ_FILE){
// 		cerr<<"ERROR: Provide a sequence to scan with the -s option"<<endl;
// 		EXIT=true;
// 	}
// 	if(PSSM && AFFINITY_SCORE ){
// 		cerr<<"ERROR: Score matrix not applicable with affinity based annotation"<<endl;
// 		EXIT=true;
// 	}
// 	if(USER_SHIFT && !AFFINITY_SCORE){
// 		cerr<<"ERROR: Shift operation not possible while using traditional annotation"<<endl;
// 		EXIT=true;
// 	}
// 	if(FACT_LAMBDA && !AFFINITY_SCORE){
// 		cerr<<"ERROR: LAMBDA paramater cannot be used with traditional score annotaition"<<endl;
// 		EXIT=true;
// 	}
// 	if(FACT_LNR0 && !AFFINITY_SCORE){
// 		cerr<<"ERROR: LNR0 parameter cannot be used with trtadtional score annotation"<<endl;
// 		EXIT=true;
// 	}
// 	if(EXTN_GC_FLAG && REGION && !AFFINITY_SCORE){
// 		cerr<<"ERROR: Cannot specify GC global externally and local GC for regions at the same time"<<endl;
// 		EXIT=true;
// 	}

	if (MASK_SCORING) {
		if (!mask_handle.compare("skip")){
			MASK_SCORING = false;
		} else if (!mask_handle.compare("worse")){
			ANNOTATION_AVERAGE = false;
		} else if (!mask_handle.compare("average")){
			ANNOTATION_AVERAGE = true;
		} else {
			cerr<<"ERROR: the value "<<mask_handle<<" is not recognized for the option --mscore."<<endl;
			EXIT=true;
		}
	}

// 	if (DISCARD_POS) {

// 		int pos_value;
// 		int unsigned next;

// 		while( (next= discard_pos_param.find_first_of(',')) != discard_pos_param.npos ) {
// 			if (next> 0)
// 				pos_value = atoi(discard_pos_param.substr(0,next).c_str());
// 			if (pos_value == 0) {
// 				cerr<<"ERROR: The position to discard is not specified correctly"<<endl;
// 				EXIT=true;
// 			}
// 			discard_pos.push_back(pos_value-1); //-1 because the internal motif position number starts at 0
// 			discard_pos_param = discard_pos_param.substr(next+1);
// 		}
// 		if(discard_pos_param.length() > 0)
// 			pos_value = atoi(discard_pos_param.substr(0,next).c_str());
// 		if (pos_value == 0) {
// 			cerr<<"ERROR: The position to discard is not specified correctly"<<endl;
// 			EXIT=true;
// 		}
// 		discard_pos.push_back(pos_value-1); //-1 because the internal motif position number starts at 0
// 	}

	if (SINGLE_STRAND){
		if (!AFFINITY_SCORE){
				cerr<<"ERROR: the --strand option is only compatible with affinity calculation"<<endl;
				EXIT=true;
			} else {
					if (!single_strand_direct.compare("F")){
							SS_FORWARD = true;
					} else if (!single_strand_direct.compare("R")){
							SS_FORWARD = false;
					} else {
						cerr<<"ERROR: the value "<< single_strand_direct <<" is not recognized for the option --strand (can only be F or R)."<<endl;
						EXIT=true;
					}
			}
	}


	/* if(abs(atof(regSize))!=atof(regSize)){
     cerr<<"ERROR: Region size can only take real values"<<endl;
     }
     if(abs(atof(regShift))!=atof(regShift)){
     cerr<<"ERROR: Region Shift can only take real values"<<endl;
     }*/
	if(EXIT){
		// error print
		// cerr<<"Invalid options. Check manual for proper usage with -h option."<<endl;
		// exit(1);

        throw invalid_argument("ERROR: Invalid argument(s) for TRAP.");
        return;
	}
}

///////////////-----------------------------------------------------/////////////////////////
// This function was the original main function and adapted for parameters needed for PROmiRNA
// If more parameters are required they need to be added in the future, the code to process them still exists
void annotate(vector<double> & affinities, char * fastaPath, char* psemPath, float cpgContent)
{
    // Set global variables for our parameters
    // NOTE: Below is a list of parameters and the global variables that need to be set for them
    // This could become important in the future if more parameters need to be added

    EXTN_SEQ_FILE = true; // For FASTA input file
    seqFile = fastaPath; // For FASTA input file

    // EXTN_GFF_FILE = true; // For GFF output file
    // op_file = gffPath; // For GFF output file

    PSEM = true; // For PSEM input file
    AFFINITY_SCORE = true; // For PSEM input file
    matFile = psemPath; // For PSEM input file

    DEF_GC = true; // For CpG content
    EXTN_GC_FLAG = true; // For CpG content
    gc_content = cpgContent; // For CpG content

    // checking for wrong inputs during entry of parameters in the begining
    validate_options();

    // initialize the input of sequence file and output of gff file
    initialize_Sequenceinput_GFFoutput();

    // reading the matrix file.. and storing all the matrices in motif_matrices list
    read_matrix_file(matFile); //see matrix.h and read_matrix.cpp files

    // if affinity based scoring
    if(AFFINITY_SCORE)
    {
        analyze_matrix_properties();
        if (MASK_SCORING)
        {
            matrices_mask_scoring();
        }
        // start the routine for reading sequences and score calculation
        calc_affinities();
    }
    else
    {
        analyze_matrix_properties();
        // call the routine for reading sequences and score calculation
        calc_hits();
    }

    motif_matrices.clear();
    pssm_motif_matrices.clear();
    psem_motif_matrices.clear();

    affinities = tmp_affinities;
}

}

///////////////-----------------------------------------------------/////////////////////////
// From the original argument parsing function this list shows which global variables need to be adapted if more parameters are added in the future
// For our first implementation we set them for the parameters we use manually here

// Option           Variables that need to be changed
// -s               EXTN_SEQ_FILE = true, seqFile = given parameter (char *)
// -o               EXTN_GFF_FILE = true, op_file = given parameter (char *)
// -R               REGION = true, regSize = given parameter (string), region_size = given parameter (int)
// -t               EXTN_THRESH = true
// -ttype           THRESH_TYPE = true
// -g               DEF_GC = true, EXTN_GC_FLAG = true, gc_content = given parameter (float)
// -a               AFFINITY_SCORE = true
// -d               DEBUG = true
// -l               LARGEST = true, nTopScore = given parameter (string), number_of_top_scores = given parameter (int)
// -S               USER_SHIFT = true, regShift = given parameter (string), region_shift = given parameter (int)
// -lambda          FACT_LAMBDA = true, fact_lambda = given parameter (float)
// -lnR0            FACT_LNR0 = true, fact_lnR0 = given parameter (float)
// -pssm            PSSM = true, matFile = given parameter (string)
// -psem            PSEM = true, AFFINITY_SCORE = true, matFile = given parameter (string)
// -mscore          mask_handle = given parameter (string)
// -discard         DISCARD_POS = true, discard_pos_param = given parameter (string)
// -tab             TABLE = true
// -strand          SINGLE_STRAND = true, single_strand_direct = given parameter (string)
// -name            NAME_FLAG = true, char * name_fl = given parameter, char * names = given parameter split into tokens (" ,"-separated) afterwards those are written into vector "name_flags"
// -gap             GAP = TRUE (cannot be more than one character)
