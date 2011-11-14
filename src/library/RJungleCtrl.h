#ifndef RANDOMJUNGLECTRL_H_
#define RANDOMJUNGLECTRL_H_

#ifdef __DEBUG__
#define RANDOMJUNGLECTRLDEBUG
#endif

/*
 * Includes
 */

#include <iostream>
#include <cstring>
#include <vector>
#include <string>
#include <limits>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "treedefs.h"

// base headers
#include "RJungleIO.h"
#include "RJunglePar.h"
#include "RJungleGen.h"
#include "RJungleFromXML.h"
#include "RJungleCompiler.h"
#include "RJungleConfusion.h"
#include "RJunglePrediction.h"
#include "RJungleImportance.h"

// enhanced headers
#include "RJungleBinder.h"
#include "RJungleGrow.h"
#include "RJungleProxi.h"
#include "RJungleImpute.h"
#include "RJungleTuneMtry.h"
#include "RJungleProto.h"
#include "RJungleOutlier.h"

#include "TimeProf.h"

/*
 * Def.
 */

template<class T>
class RJungleCtrl {
public:
  RJungleCtrl() {
  }

  virtual ~RJungleCtrl() {
  }

  void setParAndIO(RJunglePar &par_, RJungleIO &io_) {
    binder.par = par_;
    binder.io = io_;
  }

  /*
   * internal method:
   * automatically make the jungle
   */
  void autoBuildInternal(RJunglePar &par_, RJungleIO &io_, RJungleGen<T> &gen_,
      DataFrame<T> &data, std::vector<uli_t> *&colMaskVec) {

    try {
      binder.par = par_;
      binder.io = io_;
      binder.gen = gen_;

      // importance
      binder.intrinsicImportance = gen_.fct.newImportanceObject();
      binder.permImportance = new PermImportance<T> ();

      binder.intrinsicImportance->setBaseData(&par_, &io_, &data);
      binder.permImportance->setBaseData(&par_, &io_, &gen_, &data);

      std::vector<std::pair<double, uli_t> >::reverse_iterator itFreq;
      std::map<uli_t, double>::iterator itMap;
      std::map<uli_t, double>::iterator itMapD;

      DataFrame<double> *yaimpVariation = NULL;

      int k, totalTrees;
      int i = 0;
      int iterImpute = 0;
      double timeEst = 0;
      bool mtryIsUndef = (binder.par.mtry == 0) ? true : false;
      bool ntreeAutoTweek = (binder.par.ntree == 0) ? true : false;
      double probNtree = 0.999;
      double mtryRatio = (double) binder.par.mtry / binder.par.ncol;
      double divFactor = 1.0;
      uli_t numOfVars = 0;
      uli_t mtryOld = 0;

      // impute missing data
      if (binder.par.imputeIt > 0) {
        /*if (binder.par.gwa_flag)
         RJungleImpute<T >::imputeTriviallyGWA(binder.par, data);
         else
         */
        RJungleImpute<T>::imputeTrivially(data);
      }

      if (binder.par.backSel == bs_50P) {
        divFactor = 0.5;
      } else if (binder.par.backSel == bs_33P) {
        divFactor = 2.0 / 3;
      }

      if (binder.par.backSel)
        *binder.io.outVerbose << "Backward elimination:" << std::endl;

      //if (binder.par.backSel == bs_DIAZURIATE)
      //  *binder.io.outVerbose << "Diaz-uriate selection:" << std::endl;

      // backup variables
      RJunglePar parBak = binder.par;
      std::vector<uli_t> *colMaskVecBak, *colMaskVecProt;
      colMaskVecBak = colMaskVecProt = NULL;

      if (colMaskVec != NULL) {
        colMaskVecBak = new std::vector<uli_t>();
        *colMaskVecBak = *colMaskVec;
      }

      do { //imputation

        //restore variables
        if (colMaskVecBak != NULL) {
          *colMaskVec = *colMaskVecBak;
        } else {
          if (colMaskVec == NULL)
            delete colMaskVec;
          colMaskVec = NULL;
        }
        binder.par = parBak;

        binder.intrinsicImportance->init();
        if (RJungleGrow<T>::doPermImp(binder.par, colMaskVec))
          binder.permImportance->init();

        do { // jungle growing
          // reset cache
          binder.intrinsicImportance->reset();
          if (RJungleGrow<T>::doPermImp(binder.par, colMaskVec))
            binder.permImportance->reset();

          // burning down the jungle
          if (binder.trees.size() != 0) {
            typename std::vector<CmpldTree<T> *>::iterator it;
            it = binder.cmpldTrees.begin();
            while (it != binder.cmpldTrees.end()) {
              delete *it;
              *it = NULL;
              ++it;
            }
          }

          // create yaimp
          if (binder.par.yaimp_flag) {
            /*
             binder.yaimp = new DataFrame<double >(binder.par);
             binder.yaimp->setDim(binder.par.ncol, binder.par.ncol);
             binder.yaimp->initMatrix();
             binder.yaimp->setVarNames(data.varNames);
             binder.yaimp->setAll(0);

             yaimpVariation = new DataFrame<double >(binder.par);
             yaimpVariation->setDim(binder.par.ncol, binder.par.ncol);
             yaimpVariation->initMatrix();
             yaimpVariation->setVarNames(data.varNames);
             yaimpVariation->setAll(0);
             */
          }

          k = 0;

          // num of rand. selected var.
          if (mtryIsUndef) {
            binder.par.mtry = binder.gen.fct.setMtry((colMaskVec)
                ? (colMaskVec->size()) : (binder.par.ncol - 1));

          } else if (colMaskVec != NULL) {
            binder.par.mtry = (uli_t) (colMaskVec->size() * mtryRatio);
          }

          // correct mtry
          if (binder.par.mtry > binder.par.ncol - 1)
            binder.par.mtry = binder.par.ncol - 1;

          //if (binder.par.mtry < 2)
          //  binder.par.mtry = 2;

          if (binder.par.mtry < 1)
            binder.par.mtry = 1;

          // get new tree size
          if (ntreeAutoTweek) {
            if (colMaskVec == NULL) {
              // set inital treesize
              binder.par.ntree = Helper::getTreeSizeProb(probNtree,
                  binder.par.mtry, binder.par.ncol);
              this->binder.par.ntree = binder.par.ntree;
            } else {
              binder.par.ntree = Helper::getTreeSizeProb(probNtree,
                  binder.par.mtry, colMaskVec->size());
              this->binder.par.ntree = binder.par.ntree;
            }
          }
          mtryOld = binder.par.mtry;

          //if (binder.par.mtry < 2)
          //  throw Exception(ERRORCODE_33);
          assert(binder.par.mtry > 0);

          if (binder.par.ncol - 1 < binder.par.mtry)
            throw Exception(ERRORCODE_38);

          //set tree vector
          if (binder.par.ntree == 0)
            throw Exception(ERRORCODE_32);

        uli_t myTreeSize = binder.par.ntree;
#ifdef HAVE_MPI
        if (binder.par.mpiId == 0)
          myTreeSize = binder.par.ntreeMpi;
#endif
        binder.oobSet.init(binder.par.nrow, myTreeSize);

#ifdef __DEBUG2__
    //*io.outVerbose << data << std::endl;
#endif
          // grow trees
          *binder.io.outVerbose << "Number of variables: " << ((colMaskVec
              == NULL) ? (binder.par.ncol) : (colMaskVec->size()));
          *binder.io.outVerbose << " (mtry = " << binder.par.mtry << ")"
              << std::endl;

          if (binder.par.backSel) {
            *binder.io.outVerbose << "Targeted number of variables <= "
                << binder.par.numOfImpVar << std::endl << "Iteration: " << i
                << std::endl;
          }

          timeEst = 0;
          totalTrees = 0;

					// grow trees
          RJungleGrow<T>::growLocal(data, binder, colMaskVec, i);

          if (binder.par.yaimp_flag) {
            //calc. z-score
            /*
             yaimpVariation->minusXPow2(*binder.yaimp);
             yaimpVariation->div((double)binder.par.ntree);
             yaimpVariation->squareroot();
             binder.yaimp->div(*yaimpVariation);

             binder.yaimp->printCSV(*binder.io.outYaimp);
             */
          }

          // save column vector
          if ((binder.par.prototypes > 0) || (binder.par.varproximities > 0)) {
            if (colMaskVec != NULL) {
              if (colMaskVecProt == NULL)
                colMaskVecProt = new std::vector<uli_t>();

              *colMaskVecProt = *colMaskVec;
            }
          }

          // output intrinsic importance
          binder.intrinsicImportance->setIteration(i);

#ifdef __DEBUG__
        *binder.io.outVerbose << "Save gini (if mpiId >= 1)." << std::endl;
#endif
        
#ifdef HAVE_MPI
          if (binder.par.mpiId == 0)
            binder.intrinsicImportance->save();
#else
          binder.intrinsicImportance->save();
#endif

#ifdef __DEBUG__
        *binder.io.outVerbose << "Eval. perm. imp. results" << std::endl;
#endif
          // evaluate permutation importance results
          if (RJungleGrow<T>::doPermImp(binder.par, colMaskVec)) {

            // discard variables with low importance values
            if (colMaskVec == NULL) {
              numOfVars = binder.par.ncol;
              colMaskVec = new std::vector<uli_t>();
            } else {
              numOfVars = colMaskVec->size();
              colMaskVec->clear();
            }

            binder.permImportance->getBestVariables((uli_t) floor(numOfVars
                  * divFactor), *colMaskVec);
          } else {

            if (colMaskVec == NULL) {
              numOfVars = binder.par.ncol;
              colMaskVec = new std::vector<uli_t>();
            } else {
              numOfVars = colMaskVec->size();
              colMaskVec->clear();
            }

            binder.intrinsicImportance->getBestVariables((uli_t) floor(
                numOfVars * divFactor), *colMaskVec);
          }

#ifdef HAVE_MPI
#ifdef __DEBUG__
        *binder.io.outVerbose << "send/receive colMaskVec" << std::endl;
#endif
          // broadcast new colMaskVec to mpi slaves
          if (binder.par.mpi == 1) {
            uli_t iLength;

            if (binder.par.mpiId == 0) {
              // master sends data

              // send length of vector
              iLength = colMaskVec->size();
              MPI_Bcast(&iLength, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

              // send vector
              MPI_Bcast(&(colMaskVec->at(0)), iLength, MPI_UNSIGNED, 0,
                  MPI_COMM_WORLD);

              // send savecollector
              MPI_Bcast(&(binder.permImportance->saveCol.isAvailable[0]),
                  binder.par.ncol, MPI_CHAR, 0, MPI_COMM_WORLD);

            } else {
              // slaves receive data

              // receive length
              MPI_Bcast(&iLength, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
              if (colMaskVec == NULL)
                colMaskVec = new std::vector<uli_t>();
              colMaskVec->resize(iLength);

              // receive vectors
              MPI_Bcast(&(colMaskVec->at(0)), iLength, MPI_UNSIGNED, 0,
                  MPI_COMM_WORLD);

              MPI_Bcast(&(binder.permImportance->saveCol.isAvailable[0]),
                  binder.par.ncol, MPI_CHAR, 0, MPI_COMM_WORLD);
            }
          }
#endif
          //          *binder.io.outVerbose << "jo (" << std::endl;
//          for (uli_t jjj = 0; jjj < colMaskVec->size(); ++jjj)
//            *binder.io.outVerbose << (double) colMaskVec->at(jjj) << " ";
//          *binder.io.outVerbose << std::endl;

          //          *binder.io.outVerbose << ")" << std::endl;

          if ((binder.par.backSel == bs_DIAZURIATE) && (colMaskVec != NULL)) {
            numOfVars = (colMaskVec == NULL) ? binder.par.ncol
                : colMaskVec->size();
          }

          // iteration

#ifdef __DEBUG__
          *binder.io.outVerbose << "Print conf. matrix" << std::endl;
#endif
          // print confusion matrix
					TIMEPROF_START("RJungleCtrl~~RJungleConfusion::getAndPrintOobPredictions");
          RJungleConfusion<T>::printConfMatNew(binder.io, binder.gen, data,
              binder.pred.pred, binder.oobSet, i, numOfVars);
					TIMEPROF_STOP("RJungleCtrl~~RJungleConfusion::getAndPrintOobPredictions");

					// Write OOB/tree data
					if (binder.par.oob_flag)
						*binder.io.outOOB << binder.oobSet << std::endl << std::endl;

          // stop, if no back. sel.
          if (binder.par.backSel == 0)
            break;

          // do backward elimination
          if ((binder.par.backSel != bs_DIAZURIATE) || ((binder.par.backSel
              == bs_DIAZURIATE) && (colMaskVec == NULL))) {

            // stop criteria for backward elimination
            if ((binder.par.backSel != bs_DIAZURIATE) && (numOfVars
                <= binder.par.numOfImpVar))
              break;
          }

          if (binder.par.backSel == bs_DIAZURIATE) {
            // save numOfVars

            // delete to threshold binder.par.numOfImpVar
            //do {
            colMaskVec->pop_back();
            //} while (colMaskVec->size() > binder.par.numOfImpVar);

            // get accuracy of this step

            //*binder.io.outVerbose << "Num. of var.: " << colMaskVec->size()
            //  << " accuracy: " << binder.gen.fct.getAccuracyOfCmpldTrees(
            //  &data, binder.cmpldTrees, binder.gen.fct, data.getDepVar(),
            //  binder.oobSet, -1, false, binder.io);

            // stop criteria
            if (colMaskVec->size() == 4)
              break;
          }

          // clean up
          if (binder.par.yaimp_flag) {
            delete binder.yaimp;
            binder.yaimp = NULL;
            delete yaimpVariation;
            yaimpVariation = NULL;
          }

          *binder.io.outVerbose << std::endl;

          i++;
        } while (true);

        // impute data in the first step
        if (binder.par.imputeIt > 0) {
          ++iterImpute;
          if (binder.par.imputeIt == (unsigned int) iterImpute) {
            *binder.io.outVerbose << "Writing imputed data to file..."
                << std::endl;
            sort(colMaskVec->begin(), colMaskVec->end());
            data.printCSV(*binder.io.outImputedData, colMaskVec);
            break;
          } else {
            *binder.io.outVerbose << "Imputing data..." << std::endl;

            // impute data
            RJungleImpute<T>::imputeData(binder.par, binder.io, data,
                colMaskVec, binder.proxi);

            *binder.io.outVerbose << std::endl;
            continue;
          }
        } else {
          break;
        }
      } while (true);

      // get outlier scores
      if (Helper::isClassTree(binder.par) && (binder.par.outlier > 0)) {
        RJungleOutlier::getOutlier(binder.io, data, binder.proxi);
      }

      // print sample proximities
      if (binder.par.sampleproximities_flag) {
        *binder.io.outVerbose << "Writing sample proximities to file..."
            << std::endl;
        binder.proxi.printCSV(*binder.io.outSamProximity);
      }

      // print variable proximities

      // print prototypes
      if (Helper::isClassTree(binder.par) && (binder.par.prototypes > 0)) {
        *binder.io.outVerbose << "Getting prototypes..." << std::endl;
        RJungleProto<T>::getPrototypes(binder.io, binder.gen, data,
            colMaskVecProt, binder.proxi);
      }

      if (colMaskVec != NULL) {
        delete colMaskVec;
        colMaskVec = NULL;
      }

      if (colMaskVecBak != NULL) {
        delete colMaskVecBak;
        colMaskVecBak = NULL;
      }

      if (colMaskVecProt != NULL) {
        delete colMaskVecProt;
        colMaskVecProt = NULL;
      }
    } catch (std::exception &e) {
      std::string *out = new std::string("RJungleCtrl::make:");
      out->append(e.what());
      throw Exception(out->c_str());
    }
  }

  void imputeSNPs(RJunglePar &par, RJungleIO &io, RJungleGen<T> &gen,
      DataFrame<T> &data, std::vector<uli_t> *colMaskVec) {

    RJungleImpute<T>::imputeSNPs(par, io, gen, data, colMaskVec);
  }

  void tuneMtry(RJunglePar & par, RJungleIO & io, RJungleGen<T> & gen,
      DataFrame<T> & data, std::vector<uli_t> *colMaskVec) {

    RJungleTuneMtry<T>::tuneMtry(par, io, gen, data, colMaskVec);
  }

  /*
   * automatically make the jungle with parameters given in "par"
   */
  static void autoBuild(RJunglePar par) {
    std::vector<uli_t> *colMaskVec = NULL;
    time_t start, end;
    clock_t startgrow, endgrow;

    // validate par
    if ((par.backSel > 0) && (par.tunemtry != 0))
      throw Exception(ERRORCODE_51);

    if ((par.backSel > 0) && (par.sampleproximities_flag != 0))
      throw Exception(ERRORCODE_52);

    if ((par.backSel > 0) && (par.varproximities != 0))
      throw Exception(ERRORCODE_53);

    // create controller
    RJungleCtrl<T> rjCtrl;

    // input output
    if (strcmp(par.outprefix, "") == 0)
      par.outprefix = (char *) "rjungle";
    RJungleIO io;
    io.open(par);

    // header
    time(&start);
    RJungleHelper<T>::printHeader(par, io, start);

    /* init parameters:
     * M. Matsumoto and T. Nishimura, ACM Transactions on Modeling and
     * Computer Simulation, Vol. 8, No. 1 (Jan. 1998), Pages 3-30
     */
    par.rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(par.rng, par.seed);
    Helper::getCSVdim(par);
    *io.outVerbose << "Read " << par.nrow << " row(s) and " << par.ncol
        << " column(s)." << std::endl;

    // init data frame
    DataFrame<T> *data = new DataFrame<T> (par);
    RJungleGen<T> gen;

    // get data from file
    par = data->getDataFromCSV(); // and return updated parameters
#ifdef __DEBUG__
    //*io.outVerbose << *data << std::endl;
#endif
		
    // read parameter if prediction is wanted
		if (strcmp(par.predict, "") != 0) {
			RJungleFromXML<T> xmlReader(par);
			xmlReader.data = data;
			par = xmlReader.getParFromXmlFile();
		}

    // permute response?
    if (par.permresponse_flag) {
      T *vec = new T[par.nrow];
      T *bakVec = new T[par.nrow];

      *io.outVerbose << "Permute response..." << std::endl;
      data->getRowPieceWise(par.depVar, vec);
      gsl_ran_shuffle(par.rng, vec, par.nrow, sizeof(T));
      data->setCol(par.depVar, vec);

      delete[] vec;
      delete[] bakVec;
    }

    *io.outVerbose << "Use " << par.nrow << " row(s) and " << par.ncol
        << " column(s)." << std::endl;

    // GeneratorFct
    gen.init(par, *data);

    // get column selection
    if (strcmp(par.colSelection, "") != 0) {
      *io.outVerbose << "Getting selected columns (variables) from "
          << par.colSelection << "..." << std::endl;

      std::vector<std::string> selNames;
      selNames = Helper::getColSelection(par.colSelection, par.delimiter);
      colMaskVec = data->getIndexOfVarNames(selNames);
    }

    // ...
    if (par.extractdata_flag) {
      *io.outVerbose << "Extracting data..." << std::endl;
      data->printCSV(*io.outExtractData, colMaskVec);
    } else if (par.tunemtry > 0.0) {
      *io.outVerbose << "Tuning mtry..." << std::endl;
      rjCtrl.tuneMtry(par, io, gen, *data, colMaskVec);
    } else {

      *io.outVerbose << "Dependent variable name: " << par.depVarName
          << std::endl;
      /*
       if (par.gwa_flag && par.imputeIt > 0) {
       // impute SNP data
       *io.outVerbose << "Imputing SNP data..." << std::endl;
       rjCtrl.imputeSNPs(par, io, gen, *data, colMaskVec);
       } else {
       // automatically grow jungle
       *io.outVerbose << "Growing jungle..." << std::endl;
       rjCtrl.autoBuildInternal(par, io, gen, *data, colMaskVec);
       }
       */

      *io.outVerbose << "Growing jungle..." << std::endl;
#ifdef HAVE_MPI
      *io.outVerbose << "Using " << par.mpiSize << " process(es)." << std::endl;
#endif

      startgrow = clock();
			TIMEPROF_START("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
      rjCtrl.autoBuildInternal(par, io, gen, *data, colMaskVec);
			TIMEPROF_STOP("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
      endgrow = clock();
    }

    // print info stuff
    RJungleHelper<T>::printRJunglePar(par, *io.outLog);

    // clean up
    if (data != NULL)
      delete data;
    if (colMaskVec != NULL)
      delete colMaskVec;

    time(&end);
    RJungleHelper<T>::printFooter(par, io, start, end, startgrow, endgrow);

    // free random number generator
    gsl_rng_free(par.rng);

#ifdef HAVE_TIMEPROF
		timeProf.writeToFile(*io.outTimeProf);
#endif

    io.close();
  }

private:
  // binder contains all important variables
  RJungleBinder<T> binder;
};

#endif /*RANDOMJUNGLECTRL_H_*/
