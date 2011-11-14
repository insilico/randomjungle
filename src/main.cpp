/* This software is called Random Jungle (RJ) and is a fast implementation
 * of Random Forests. The software implements several features
 * such as variable importance or classifier creation.
 * Find software at http://www.randomjungle.org webpage.
 * Copyright (C) 2008-2011  Daniel F. Schwarz
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"
#include <iostream>
#include <exception>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include "lr.h"
#include "librjungle.h"

/*
 * Source
 */

int flag[1];

/*! 
 * \brief Main function handles arguments and calls starting point in rjungle lib.
 * 
 * @param argc Number of arguments
 * @param argv Arguments
 * @return Error code
 */
int main(int argc, char **argv) {
  RJunglePar par;
	int c;

  par = initRJunglePar();

	while (1) {
    //abcdefghijklmnopqrstuvwxyz
    //ABCDEfGHI KLMNOPQRSTUVWXYZ
	  char const *argstring =
	    "f:e:y:t:m:d:r:c:n:g:o:a:b:x:i:j:hvuw:pVSP:C:D:I:M:GATsWEYz:U:X:Zk:l:B:Q:q:O:N:H:K:L";

    #ifdef HAVE_GETOPT_LONG
	  static struct option long_options[] = {
      {"file",                required_argument, 0, 'f'},
      {"delimiter",           required_argument, 0, 'e'},
      {"treetype",            required_argument, 0, 'y'},
      {"ntree",               required_argument, 0, 't'},
      {"mtry",                required_argument, 0, 'm'},
      {"depvar",              required_argument, 0, 'd'},
      {"depvarname",          required_argument, 0, 'D'},
      {"nrow",                required_argument, 0, 'r'},
      {"ncol",                required_argument, 0, 'c'},
      {"varnamesrow",         required_argument, 0, 'n'},
      {"depvarcol",           required_argument, 0, 'g'},
      {"outprefix",           required_argument, 0, 'o'},
      {"skiprow",             required_argument, 0, 'a'},
      {"skipcol",             required_argument, 0, 'b'},
      {"missingcode",         required_argument, 0, 'x'},
      {"impmeasure",          required_argument, 0, 'i'},
      {"backsel",             required_argument, 0, 'B'},
      {"nimpvar",             required_argument, 0, 'j'},
      {"help",                no_argument, 0, 'h'},
      {"downsampling",        no_argument, 0, 'u'},
      {"verbose",             no_argument, 0, 'v'},
      {"memmode",             required_argument, 0, 'M'},
      {"write",               required_argument, 0, 'w'},
      {"predict",             required_argument, 0, 'P'},
      {"varproximities",      required_argument, 0, 'H'},
      {"summary",             no_argument, 0, 'S'},
      {"testlib",             no_argument, &flag[0], 1},
      {"plugin",              required_argument, 0, 'X'},
      {"colselection",        required_argument, 0, 'C'},
      {"impute",              required_argument, 0, 'I'},
      {"gwa",	                no_argument, 0, 'G'},
      {"impcont",	            no_argument, 0, 'A'},
      {"transpose",	          no_argument, 0, 'T'},
      {"sampleproximities",   no_argument, 0, 's'},
      {"weightsim",           no_argument, 0, 'W'},
      {"extractdata",         no_argument, 0, 'E'},
      {"yaimp",               no_argument, 0, 'Y'},
      {"seeed",               required_argument, 0, 'z'},
      {"nthreads",            required_argument, 0, 'U'},
      {"pedfile",             no_argument, 0, 'p'},
      {"version",             0, 0, 'Z'},
      {"maxtreedepth",        required_argument, 0, 'k'},
      {"targetpartitionsize", required_argument, 0, 'l'},
      {"pluginpar",           required_argument, 0, 'Q'},
      {"tunemtry",            required_argument, 0, 'q'},
      {"outlier",             required_argument, 0, 'O'},
      {"votes",               no_argument, 0, 'V'},
      {"prototypes",          required_argument, 0, 'N'},
      {"condimp",             required_argument, 0, 'K'},
      {"permresponse",        no_argument, 0, 'L'},
      {"oobset",              no_argument, 0, 1000},
      {"classweights",        required_argument, 0, 1001},
      {0, 0, 0, 0}};
    int option_index = 0;

    c = getopt_long(
      argc,
      argv,
      argstring,
      long_options,
      &option_index);

    #else
    c = getopt(
      argc,
      argv,
      argstring
      );
    #endif

    // Detect the end of the options.
    if (c == -1) break;

    switch (c) {
      case 0:
        // If this option set a flag, do nothing else now.
        #ifdef HAVE_GETOPT_LONG
        if (long_options[option_index].flag != 0) break;
        std::cout << "option " << long_options[option_index].name;
        #endif
        if (optarg) std::cout << " with arg " << optarg;
        std::cout << std::endl;
        break;

  		case 'f':
  			par.filename = optarg;
  			break;

  		case 'e':
  			par.delimiter = optarg[0];
  			break;

  		case 'y':
  			par.treeType = atoi(optarg);
  			break;

  		case 't':
  			par.ntree = atoi(optarg);
  			break;

  		case 'm':
  			par.mtry = atoi(optarg);
  			break;

  		case 'd':
  			par.depVar = atoi(optarg);
  			break;

  		case 'r':
  			par.nrow = atoi(optarg);
  			break;

  		case 'c':
  			par.ncol = atoi(optarg);
  			break;

  		case 'n':
  			par.varNamesRow = atoi(optarg);
  			break;

  		case 'g':
  			par.depVarCol = atoi(optarg);
  			break;

  		case 'a':
  			par.skipRow= atoi(optarg);
  			break;

  		case 'b':
  			par.skipCol = atoi(optarg);
  			break;

  		case 'x':
  			par.missingcode = atoi(optarg);
  			break;

  		case 'o':
  			par.outprefix = optarg;
  			break;

      case 'i':
        par.impMeasure = atoi(optarg);
        break;

      case 'B':
        par.backSel = atoi(optarg);
        break;

  		case 'j':
  			par.numOfImpVar = atoi(optarg);
  			break;

      case 'v':
  			par.verbose_flag = 1;
        break;

      case 'u':
  			par.downsampling_flag = 1;
        break;

      case 'w':
  			par.saveJungleType = atoi(optarg);
        break;

      case 'P':
  			par.predict = optarg;
        break;

      case 'H':
  			par.varproximities = atoi(optarg);
        break;

      case 'S':
        par.summary_flag = 1;
        break;

      case 'M':
  			par.memMode = atoi(optarg);
      break;

      case 'X':
  			par.plugin = optarg;
        break;

      case 'C':
  			par.colSelection = optarg;
        break;

      case 'D':
  			par.depVarName = optarg;
        break;

      case 'I':
        par.imputeIt = atoi(optarg);
        break;

      case 'G':
        par.gwa_flag = 1;
        break;

      case 'A':
        par.allcont_flag = 1;
        break;

      case 'T':
        par.transpose_flag = 1;
        break;

      case 's':
        par.sampleproximities_flag = 1;
        break;

      case 'W':
        par.weightsim_flag = 1;
        break;

      case 'E':
        par.extractdata_flag = 1;
        break;

      case 'Y':
        par.yaimp_flag = 1;
        break;

      case 'z':
        par.seed = atoi(optarg);
        break;

      case 'U':
        par.nthreads = atoi(optarg);
        break;

      case 'p':
        par.pedfile_flag = 1;
        break;

      case 'Z':
        std::cout << "Random Jungle (RJ) version: " << par.version <<
#ifdef HAVE_MPI
					" (mpi)" <<
#endif		
				  std::endl;
        return 0;
        break;

      case 'k':
        par.maxTreeDepth = atoi(optarg);
        break;

      case 'Q':
        par.pluginPar = optarg;
        break;

      case 'l':
        par.targetPartitionSize = atoi(optarg);
        break;

      case 'q':
        par.tunemtry = atof(optarg);
        break;

      case 'K':
        par.condimp = atof(optarg);
        break;

      case 'O':
        par.outlier = atoi(optarg);
        break;

      case 'V':
        par.votes_flag = 1;
        break;

      case 'N':
        par.prototypes = atoi(optarg);
        break;

      case 'L':
        par.permresponse_flag = 1;
        break;

			case 1000:
				par.oob_flag = 1;
				break;

			case 1001:
				par.classweights = optarg;
				break;

      case 'h':
        std::cout << "Usage: jungle.exe [OPTION]..." << std::endl
        << "Analysis of data with method RandomJungle." << std::endl
        << std::endl
        << "  -f, --file=NAME         NAME of file with data (CSV/can be gzipped)" << std::endl
        << "                          EXAMPLE: data.csv, mydata.csv.gz" << std::endl
        << "                          DEFAULT: data.csv.gz" << std::endl
        << "  -o, --outprefix=NAME    NAME of output files prefixes" << std::endl
        << "                          DEFAULT: jungle" << std::endl
        << "  -e, --delimiter=CHAR    delimiter in data file is CHAR." << std::endl
        << "                          DEFAULT: ;" << std::endl
        << "  -y, --treetype=ID       choose base classifier ID:" << std::endl
        << "                          ID = 1: CART, y nominal, x numeric" << std::endl
        << "                                  CART might give a biased variable selection" << std::endl
        << "                                  on data sets with varying category sizes" << std::endl
        << "                                  per variables." << std::endl
        << "                          ID = 2: CART, y nominal, x nominal" << std::endl
        << "                                  Twoing trees." << std::endl
        << "                          ID = 3: CART, y numeric, x numeric" << std::endl
        << "                                  Regression trees." << std::endl
        << "                          ID = 4: CART, y numeric, x nominal" << std::endl
        << "                                  Regression trees." << std::endl
        << "                          ID = 5: CART2, y nominal, x numeric" << std::endl
        << "                                  CART2 is slower than CART (ID=1)" << std::endl
        << "                                  if category size is small or it is not a numeric" << std::endl
        << "                                  variable. Original Breiman/Cutler/Friedman idea." << std::endl
//
//        << "                          ID = 5: LOTUS, ordered input, categorical output," << std::endl
//        << "                                         unbiased variable selection." << std::endl
//        << "                                  LOTUS might bring small trees on noisy" << std::endl
//        << "                                  data. (unstable in program yet)" << std::endl
//        << "                          ID = 6: SCHWARZCART, ordered input, categorical output" << std::endl
//        << "                                  SCHWARZCART serve a less biased variable selection" << std::endl
//        << "                                  with permutation importance." << std::endl
//        << "                          ID = 7: MULTICART, ordered input, categorical output" << std::endl
//        << "                                  A Node may have more then two children." << std::endl
//        << "                                  Splitt points equals category size - 1." << std::endl
//        << "                          ID = 101: Conditional Inference Trees for y (nominal), x (numeric)" << std::endl
//        << "                                   (T. Hothorn et al. Unbiased Recursive Partitioning 2006 JCGS)" << std::endl
//        << "                          ID = 111: Conditional Inference Trees for SNPs with additive model." << std::endl
//        << "                                   y (nominal && two-sample), x (numeric)" << std::endl
//
        << "                          DEFAULT: 1" << std::endl

        << "  -t, --ntree=SIZE        number (SIZE) of trees in jungle. SIZE=[1-...]" << std::endl
        << "                          If SIZE==0 then the size will be set automatically" << std::endl
        << "                          depending on mtry and variable size." << std::endl
        << "                          DEFAULT: 500" << std::endl
        << "  -m, --mtry=SIZE         size (SIZE) of randomly choosen variable sets." << std::endl
        << "                          At each node building step, a variable will be" << std::endl
        << "                          selected out of the set, that serves the" << std::endl
        << "                          biggest information gain." << std::endl
        << "                          The bigger SIZE is set, the higher computing" << std::endl
        << "                          time might be." << std::endl
        << "                          The bigger SIZE is set, the more similar trees" << std::endl
        << "                          in jungle will be. SIZE=[1-...]" << std::endl
        << "                          DEFAULT: sqrt(ncol)" << std::endl
        << "  -d, --depvar=POS        output variable at column POS in the data SET!" << std::endl
        << "                          POS=[0-...]" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -g, --depvarcol=POS     output variable is the variable POS in the data FILE!" << std::endl
        << "                          POS=[0-...]" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -n, --varnamesrow=POS   variable names at row POS in data set. POS=[0-...]" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -r, --nrow=SIZE         number (SIZE) of samples. SIZE=[1-...]" << std::endl
        << "                          SIZE = 0: use all samples in data file" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -c, --ncol=SIZE         number (SIZE) of input variables (features)." << std::endl
        << "                          SIZE=[1-...]" << std::endl
        << "                          SIZE = 0: use all variables in data file" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -a, --skiprow=SIZE      skip SIZE rows (samples) before reading" << std::endl
        << "                          data from file. SIZE=[1-...]" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -b, --skipcol=SIZE      skip SIZE columns (variables) before reading" << std::endl
        << "                          data from file. SIZE=[1-...]" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -x, --missingcode=NUM   missings are coded as NUM in data file/set." << std::endl
        << "                          Advice: SIZE=[-127..127]" << std::endl
        << "                          DEFAULT: -99" << std::endl
        << "  -i, --impmeasure=ID     choose importance method." << std::endl
        << "                          ID = 1: Intrinsic Importrance (i.e. Gini index)" << std::endl
        << "                          ID = 2: Permutation Importance Breiman, Cutler (Fortran)" << std::endl
        << "                          ID = 3: Permutation Importance Liaw, Wiener (R)" << std::endl
        << "                          ID = 4: Permutation Importance, raw values, no normalization" << std::endl
        << "                          ID = 5: Permutation Importance Meng" << std::endl
        << "                          DEFAULT: 1" << std::endl
        << "  -K, --condimp=NUM       Perform conditional importance if option -i > 1." << std::endl
        << "                          NUM is the pearson's cor. coef. cutoff. The smaller NUM, the bigger" << std::endl
        << "                          a conditional importance permutation group will be created." << std::endl
        << "                          (=> More accurate, but slower)" << std::endl
        << "                          Requires: 0 <= NUM <= 1" << std::endl
        << "                          NUM < 0 => switched off" << std::endl
        << "                          DEFAULT: switched off" << std::endl
        << "  -B, --backsel=ID        choose backward elimination." << std::endl
        << "                          ID = 0: No backward elimination" << std::endl
        << "                          ID = 1: Backward elimination. Discard 50% at each step" << std::endl
        << "                          ID = 2: Backward elimination. Discard 33% at each step" << std::endl
        << "                          ID = 3: Backward elimination." << std::endl
        << "                                  Discard only negative values at each step" << std::endl
        << "                          ID = 4: D�az-Uriarte variable selection." << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -j, --nimpvar=SIZE      only necessary if backsel>0. SIZE=[1-...]" << std::endl
        << "                          how many variable should remain." << std::endl
        << "                          The lesser SIZE is, the reliable the result might be." << std::endl
        << "                          The smaller SIZE is, the higher computing time will be." << std::endl
        << "                          DEFAULT: 1" << std::endl
        << "  -v, --verbose           print some information to STDOUT" << std::endl
        << "                          DEFAULT: switched off" << std::endl
        << "  -u, --downsampling      Choose randomly samples without replacement." << std::endl
        << "                          DEFAULT: switched off (Bootstrap)" << std::endl;

        #ifndef __SPARSE_DATA__
        std::cout
        << "  -M, --memmode=ID        Usage of the heap memory (RAM)." << std::endl
        << "                          ID = 0: double precision floating point (BIG)" << std::endl
        << "                          ID = 1: single precision floating point (Normal)" << std::endl
        << "                          ID = 2: char (small)" << std::endl
        << "                                  CHAR fits normally in one byte." << std::endl
        << "                                  DATA CELL VALUE HAS TO BE" << std::endl
        << "                                  AN INTEGER IN [-127..127]" << std::endl
        << "                          DEFAULT: 0" << std::endl;
        #endif

        std::cout
        << "  -w, --write=ID          Save Random Jungle ..." << std::endl
        << "                          ID = 0: not" << std::endl
        << "                          ID = 1: to a gzipped XML file" << std::endl
        << "                          ID = 2: to a XML raw file" << std::endl
        << "                          <outprefix>.jungle.xml.gz" << std::endl
        << "                          DEFAULT: 0" << std::endl
        << "  -P, --predict=FILE      Read Random Jungle from XML file" << std::endl
        << "                          <FILE>" << std::endl
        << "                          and predict data from csv data file." << std::endl
        << "                          (parameter -f)" << std::endl
        << "                          DEFAULT: switched off (Empty String)" << std::endl
        //<< "  -H, --varproximities=ID Calculate variable proximities" << std::endl
        //<< "                          ID = 0: not" << std::endl
        //<< "                          ID = 1: Simple" << std::endl
        //<< "                          ID = 2: Extended" << std::endl
        //<< "                          DEFAULT: switched off (ID=0)" << std::endl
        << "  -S, --summary           Print a summary of all trees in jungle." << std::endl
        << "                          DEFAULT: switched off" << std::endl
        << "  -X, --plugin=FILE       Use Random Jungle generator from given plugin:" << std::endl
        << "                          <FILE>" << std::endl
        << "                          DEFAULT: switched off (Empty String)" << std::endl
        << "  -Q, --pluginpar=STRING  Additional plugin parameters." << std::endl
        << "                          DEFAULT: switched off (Empty String)" << std::endl
        << "  -C, --colselection=FILE Only use selected columns listed in <FILE>." << std::endl
        << "                          DEFAULT: switched off (Empty String)" << std::endl
        << "  -D, --depvarname=NAME   Output variable name in the data SET!" << std::endl
        << "                          DEFAULT: switched off (Empty String)" << std::endl
        << "  -I, --impute=NUM        Impute data with NUM iterations." << std::endl
        << "                          DEFAULT: switched off" << std::endl
        << "  -A, --impcont           Impute data is continuous." << std::endl
        << "                          DEFAULT: switched off (categorical data)" << std::endl

        << "  -G, --gwa               Impute data with GWA settings." << std::endl
        << "                          DEFAULT: switched off" << std::endl

        << "  -q, --tunemtry=THETA    Tune mtry parameter." << std::endl
        << "                          Step size = THETA." << std::endl
        << "                          E.g.: Quarter stepping (-q 0.25): 100 predictor variables (ncol = 100)" << std::endl
        << "                             => mtry = 25, 50, 75" << std::endl
        << "                          DEFAULT: 0 (switched off)" << std::endl
        << "  -s, --sampleproximities Calculate samples proximities." << std::endl
        << "                          DEFAULT: switched off" << std::endl
        << "  -E, --extractdata       Extract data to file." << std::endl
        << "                          DEFAULT: switched off" << std::endl
        << "  -z, --seeed=NUM         Seed of run." << std::endl
        << "                          DEFAULT: 1" << std::endl
        << "  -U, --nthreads=NUM      Use maximaly NUM threads (CPUs) for parallel processing." << std::endl
        << "                          DEFAULT: NUM OF CPUS IN COMPUTER" << std::endl
        << "  -p, --pedfile           Input file has got ped format (i.e. plink output with recodeA)." << std::endl
        << "                          DEFAULT: switched off" << std::endl
        << "  -h, --help              Displays this help text" << std::endl
        << "  -Z, --version           Displays the version number" << std::endl
        << "  -k, --maxtreedepth=NUM  Max tree depth" << std::endl
        << "                          DEFAULT: 15000" << std::endl
        << "  -l, --targetpartitionsize=NUM  Target partition size" << std::endl
        << "                          DEFAULT: 1" << std::endl
        << "  -O, --outlier=ID        Calculated outlier measure:" << std::endl
        << "                          ID = 1: raw scores Breiman's measure" << std::endl
        << "                          ID = 2: normalized Liaw's measure" << std::endl
        << "                          DEFAULT: 0 (off)" << std::endl        << std::endl
        << "  -N, --prototypes=NUM    Calculated prototypes using <NUM> nearest neighbours." << std::endl
        << "                          DEFAULT: 0 (off)" << std::endl        << std::endl
        << "  -V, --votes             Print votes for each sample." << std::endl
        << "                          DEFAULT: switched off" << std::endl
				<< "      --oobset            Output the oobset data" << std::endl
				<< "                          DEFAULT: switched off" << std::endl
				<< "      --classweights      Reweighting classes of dataset"
			                           << " (for classification trees only)" << std::endl
			  << "                          E.G.: \"1.23;22;\"" << std::endl
			  << "                                ... first class gets weight 1.23" << std::endl
			  << "                                ... second class gets weight 22" << std::endl
				<< "                          DEFAULT: switched off" << std::endl;


        #ifndef HAVE_GETOPT_LONG
        std::cout
        << "You are using a platform without function getopt_long. " << std::endl
        << "You must not use long option names (i.e. --verbose). " << std::endl
        << "Instead, use the short version (i.e. -v)." << std::endl;
        #endif

        #ifdef __SPARSE_DATA__
        std::cout
        << "This program was compiled to analyse sparse data" << std::endl
        << "with very small memory consumption." << std::endl
        << "Data values are only 0, 1, 2 or 3. (i.e. SNP data)" << std::endl
        << "If you want to use real and integer data, " << std::endl
        << "compile it without definition __SPARSE_DATA__." << std::endl;
        #endif

        return 0;
        break;

      case '?':
        //getopt_long already printed an error message.
        break;

      default:
        abort ();
    }
  }

  // Print any remaining command line arguments (not options).
  if (optind < argc) {
    printf ("non-option ARGV-elements: ");
    while (optind < argc) printf ("%s ", argv[optind++]);
    putchar ('\n');
  }

  #ifdef HAVE_MPI
  par.mpi = 1;
  MPI_Init(&argc, &argv); // all MPI programs start with MPI_Init; all 'N' processes exist thereafter
  #endif

	try {

    if (flag[0]) {
			testLibrandomjungle();
		} else {
      randomJungle(par);
		}

    #ifdef HAVE_MPI
    // mpi
    MPI_Finalize(); //  MPI Programs end with MPI Finalize; this is a weak synchronization point
    #endif

		return 0;
	} catch (std::exception &e) {

    #ifdef HAVE_MPI
    // mpi
    MPI_Finalize(); // MPI Programs end with MPI Finalize; this is a weak synchronization point
    #endif

		std::cerr 	<< "Runtime ERROR:" << e.what() << std::endl;
		return 1;
	}

}
