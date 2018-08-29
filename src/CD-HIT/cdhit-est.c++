// =============================================================================
// CD-HI-EST
// http://cd-hit.org/
// Cluster Database at High Identity (EST version)
// modified from CD-HI
//
// program written by
//                    Weizhong Li
//                    UCSD, San Diego Supercomputer Center
//                    La Jolla, CA, 92093
//                    Email liwz@sdsc.edu
//                 at
//                    Adam Godzik's lab
//                    The Burnham Institute
//                    La Jolla, CA, 92037
//                    Email adam@burnham-inst.org
//
// Modified by:
//                    Limin Fu
//                    Center for Research in Biological Systems (CRBS), UCSD
//                    La Jolla, CA, 92093
//                    Email: l2fu@ucsd.edu, fu@daovm.net
// =============================================================================

#include "cdhit-common.h"
#include "cdhit-utility.h"

namespace cdhit
{

    //over-write some defs in cd-hi.h for est version
#undef MAX_UAA
#define MAX_UAA 4

    //over-write some defs in cd-hi-init.h for est version

    void setaa_to_na();
    void make_comp_short_word_index(int NAA, int *NAAN_array, Vector<int> & Comp_AAN_idx);
    void make_comp_iseq(int len, char *iseq_comp, char *iseq);
    int cdhitEst() ;


    Options options;
    SequenceDB seq_db;

    ////////////////////////////////////  MAIN /////////////////////////////////////
    int cdhitEst(std::vector<std::string> lines, const char * outfile, unsigned const threads)
    {

        int argc = 15;
        char** argv;
        argv = (char**)calloc(argc, sizeof(char *));

        for(int i=0;i<argc;++i){
            argv[i] = (char*)calloc(256, sizeof(char));
        }
        /* THESE ABOVE need to be free'd!!!!*/
        strcpy(argv[0], "cd-hit-est");
        strcpy(argv[1] , "-o");
        strcpy(argv[2] , outfile);
        strcpy(argv[3] , "-c");
        strcpy(argv[4] , "0.8");
        strcpy(argv[5] , "-n");
        strcpy(argv[6] , "4");
        strcpy(argv[7] , "-d");
        strcpy(argv[8] , "40");
        strcpy(argv[9] , "-T");
        strcpy(argv[10] , std::to_string(threads).c_str());
        strcpy(argv[11] , "-M");
        strcpy(argv[12] , "0");

        string db_in;
        string db_out;
        string db_in_pe;
        string db_out_pe;

        //	options.cluster_thd = 0.95;
        options.NAA = 10;
        options.NAAN = NAA8;
        seq_db.NAAN = NAA8;
        options.NAA_top_limit = 12;
        setaa_to_na();
        mat.set_to_na(); //mat.set_gap(-6,-1);

        float begin_time = current_time();
        float end_time;

        // ***********************************    parse command line and open file
        if (argc < 5) print_usage_est(argv[0]);
        if (options.SetOptions( argc, argv, false, true ) == 0){
            print_usage_est(argv[0]);
        }
        options.Validate();


        for(int i=0;i<argc;++i){
            free(argv[i]);
        }
        free(argv);

        db_in     = options.input;
        db_in_pe  = options.input_pe;
        db_out    = options.output;
        db_out_pe = options.output_pe;

        InitNAA( MAX_UAA );
        seq_db.NAAN = NAAN_array[options.NAA];

        if ( options.option_r ) {
            Comp_AAN_idx.resize( seq_db.NAAN );
            make_comp_short_word_index(options.NAA, NAAN_array, Comp_AAN_idx);
        }

        if ( options.PE_mode ) {seq_db.Read( db_in.c_str(), db_in_pe.c_str(), options );}
        //else                   {seq_db.Read( db_in.c_str(),                   options, lines );}
        else                   {seq_db.DRead( db_in.c_str(),                   options, lines);}

        cout << "total seq: " << seq_db.sequences.size() << endl;
        //cout << "total seq: " << lines.size() << endl;

        seq_db.SortDivide( options );
        if(threads>1){
            seq_db.DoClustering(threads, options );
        }else{
            seq_db.DoClustering(options );
        }

        printf( "writing new database\n" );
        //if ( options.PE_mode ) {
        //seq_db.WriteClusters( db_in.c_str(), db_in_pe.c_str(), db_out.c_str(), db_out_pe.c_str(), options ); }
        // else                   {
        seq_db.WriteClusters2( db_in.c_str(),                   db_out.c_str(),                    options, lines );
        // }

        // write a backup clstr file in case next step crashes
        seq_db.WriteExtra1D( options );
        cout << "program completed !" << endl << endl;
        end_time = current_time();
        printf( "Total CPU time %.2f\n", end_time - begin_time );
        return 0;
    }

}
