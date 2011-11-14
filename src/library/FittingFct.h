#ifndef FITTINGFCT_H_
#define FITTINGFCT_H_

#ifdef __DEBUG__
//#define FITTINGFCTDEBUG
//#define FITTINGFCTDEBUG2
#endif

/*
 * Includes
 */
#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <boost/dynamic_bitset.hpp>

#include "DataFrame.h"
#include "Logistic.h"
#include "ClassAtom.h"
#include "SClassAtom.h"
#include "T2ClassAtom.h"
#include "TMClassAtom.h"
#include "IAMClassAtom.h"
#include "TermClassAtom.h"
#include "TImportance.h"
#include "IAM2WayImportance.h"
#include "Tree.h"
#include "Helper.h"
#include "INode.h"

#include "TimeProf.h"
#include "lr.h"

/*
 * Source (Def. & Decl.)
 */


template<class T>
class FittingFctPar {
public:
  FittingFctPar() :
    parent(NULL), data(NULL), depVar(0), rowMaskVec(NULL), colMaskVec(NULL),
    rowWeight(NULL), nodePerformance(0), iteration(0), depth(0), rng(NULL),
    treeImprovement(NULL), bestVar(NULL), stopGrowing(false),
        twoWayInteraction(NULL), io(NULL) {
  }

  FittingFctPar(
      INode<T> *parent, DataFrame<T> *data,
      uli_t depVar,
      std::vector<uli_t> *rowMaskVec,
      std::vector<uli_t> *colMaskVec,
      std::vector<double> *rowWeight,
      double nodePerformance,
      uli_t iteration,
      uli_t depth,
      gsl_rng *rng,
      //OUTPUT
      Importance<T> *importance,
      double *treeImprovement,
      uli_t *bestVar,
      bool stopGrowing,
      std::vector<std::pair<double, std::pair<uli_t, uli_t> > > *twoWayInteraction,
      RJungleIO *io) :
    parent(parent), data(data), depVar(depVar), rowMaskVec(rowMaskVec),
        colMaskVec(colMaskVec), rowWeight(rowWeight), nodePerformance(
            nodePerformance), iteration(iteration), depth(depth), rng(rng),
        importance(importance), treeImprovement(treeImprovement), bestVar(
            bestVar), stopGrowing(stopGrowing), twoWayInteraction(
            twoWayInteraction), io(io) {
  }

  virtual ~FittingFctPar() {
  }

  INode<T> *parent;
  DataFrame<T> *data;
  uli_t depVar;
  std::vector<uli_t> *rowMaskVec;
  std::vector<uli_t> *colMaskVec;
  std::vector<double> *rowWeight;
  double nodePerformance;
  uli_t iteration;
  uli_t depth;
  gsl_rng *rng;
  //OUTPUT
  Importance<T> *importance;
  double *treeImprovement;
  uli_t *bestVar;
  bool stopGrowing;
  std::vector<std::pair<double, std::pair<uli_t, uli_t> > > *twoWayInteraction;
  RJungleIO *io;

};

class FittingFct {
public:
  FittingFct();
  virtual ~FittingFct();
  /*
   * Get purest variable
   *
   */

  /*
   * CARTreg new
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTreg(FittingFctPar<T> &fFP) {

    // declaration and definition
    uli_t k, mt; // counter(s)
    uli_t nsp; // index of sample
    uli_t nbestt, nbest; // index of best split
    uli_t non = 0; // number of variables which did have only one kind of values
    std::vector<uli_t> &rowMask = *fFP.rowMaskVec; // selection of samples
    double &decsplit = *fFP.treeImprovement; // final decrease of sum of square
    DataFrame<T> &data = *fFP.data; // input data
    uli_t resp = fFP.depVar; // index of response variable
    T d; // response variable of a sample
    uli_t smpl; // index of sample
    uli_t kv; // index of variable
    uli_t &msplit = *fFP.bestVar; // index of best variable
    double ss; // sum of square and average value of node
    double critmax = 0; // maximal sum of square
    double crit = 0; // sum of square of split
    double ubestt = 0, ubest = 0; // result values of best split
    std::vector<uli_t> &colMask = *fFP.colMaskVec; // selection of columns
    std::vector<T> yl; // values of response variable
    std::vector<T> xt, ut; // values of indep. variable
    std::vector<std::pair<T, uli_t> > validx; // pair of variable value and rank
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> (); // classifier
    double suml, sumr, sumnode; // sums of new left, right and current node
    double npopl, npopr, nodecnt; // number of samples of new left, ...
    double av, critvar;
    av = ss = critvar = 0;

    if (rowMask.size() == 1) {// stop it, node size == 1
      fFP.stopGrowing = true;
      return classifier;
    }

    // get values of indep. variable
    data.getRowMaskedCol(resp, rowMask, yl);

    // get initial standard deviation and average of current node
    for (k = 0; k < rowMask.size(); ++k) {

      // get index of current sample
      smpl = rowMask[k];

      // get response value of current sample
      d = yl[k];

      // add amount to sum of square
      ss = ss + (double)k  * pow(av - d, 2) / (double)(k + 1);

      // add amount to average
      av = ((double)k * av + d) / (double)(k + 1);

      // get values of response variable
      yl[k] = data.at(smpl, resp);
    }

    // number of samples in current node
    nodecnt = rowMask.size();

    // sum of all indep. variable values in current node
    sumnode = av * nodecnt;

    // find best split/variable
    for (mt = 0; mt < colMask.size(); ++mt) {

      // get index of current sample
      kv = colMask[mt];

#ifdef FITTINGFCTDEBUG
      *fFP.io->outVerbose << "var: " << kv << std::endl;
#endif
      // reset decrease of sum of squares of current variable
      // critvar = 0;

      // get values of indep. variable
      data.getRowMaskedCol(kv, rowMask, xt);

      // paste indexes to indep. variable
      validx.clear();
      for (k = 0; k < rowMask.size(); ++k) {
        validx.push_back(std::make_pair(xt[k], k));
      }

      // sort indep. variable
      sort(validx.begin(), validx.end());

      // if there is only one kind of values in this variables then process next
      if ((validx.begin())->first == (validx.rbegin())->first) {
        ++non;

        if (non > 3 * data.par.ncol) { // too many variables were pure, so stop it
          fFP.stopGrowing = true;
          return classifier;
        }

        continue;
      }

      // init sums and number of samples
      suml = npopl = 0;
      sumr = sumnode;
      npopr = nodecnt;

      // get best split of current indep. variable
      for (nsp = 0; nsp < rowMask.size() - 1; ++nsp) {

        // get response value of current sample (nsp)
        d = yl[validx[nsp].second];

        // adjust sums and sample counts
        suml += d;
        sumr -= d;
        ++npopl;
        --npopr;

        if (validx[nsp].first < validx[nsp + 1].first) // value did change
          crit = (pow(suml, 2) / npopl) + (pow(sumr, 2) / npopr);

        if (crit > critvar) {// crit is the best split result til now

          // get split value
          ubestt = (validx[nsp].first + validx[nsp + 1].first) / 2.0;

          // store best value
          critvar = crit;

          // store index of sample
          nbestt = nsp;
        }
      }

      // current indep. variable had got the best split
      if (critvar > critmax) {

        // store best split value
        ubest = ubestt;

        // store best index of sample
        nbest = nbestt;

        // store index of variable
        msplit = kv;

        // store best crit value
        critmax = critvar;
      }
    }

    // store final decrease of sum of squares
    decsplit = critmax - (pow(sumnode, 2) / nodecnt);

    // build classifier
    classifier->setVarID(msplit);
    classifier->setThreshold((T) ubest);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fFP.importance)->add(msplit, decsplit / data.par.ntree);

    return classifier;
  }

  /*
   * CARTregTwoing
   *
   * Find the best nominal variable with numeric outcome.
   * The twoing method is described in CART (Leo Breiman, et al. 1983).
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTregTwoing(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    //OUTPUT
    double &decOfSos = *fitFctPar.treeImprovement;
    uli_t &sosBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    uli_t idxCell;
    T valCell;
    std::vector<T> colVec;
    std::vector<uli_t> setCur;
    std::vector<uli_t> bestSubSet;
    std::vector<uli_t> classVarVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> tmpclass;
    std::vector<uli_t> tclasspop;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxRowVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    SClassAtom<T> *classifier = new SClassAtom<T> ();
    boost::dynamic_bitset<> powerSet;

    // declaration and definition
    uli_t non = 0; // number of variables which did have only one kind of values
    std::vector<uli_t> &rowMask = *fitFctPar.rowMaskVec; // selection of samples
    uli_t resp = fitFctPar.depVar; // index of response variable
    T d; // response variable of a sample
    uli_t smpl; // index of sample
    double ss; // sum of square and average value of node
    double critmax = 0; // maximal sum of square
    std::vector<T> yl; // values of response variable
    std::vector<T> xt, ut; // values of indep. variable
    std::vector<std::pair<T, uli_t> > validx; // pair of variable value and rank
    double suml, sumr, sumnode; // sums of new left, right and current node
    double npopl, npopr, nodecnt; // number of samples of new left, ...
    double av, critvar;
    av = ss = critvar = 0;

    T val;

    if (rowMaskVec.size() == 1) {// stop it, node size == 1
      fitFctPar.stopGrowing = true;
      return classifier;
    }

    // get values of indep. variable
    data.getRowMaskedCol(resp, rowMask, yl);

    // get initial standard deviation and average of current node
    for (k = 0; k < rowMaskVec.size(); ++k) {

      // get index of current sample
      smpl = rowMaskVec[k];

      // get response value of current sample
      d = yl[k];

      // add amount to sum of square
      ss = ss + (double)k  * pow(av - d, 2) / (double)(k + 1);

      // add amount to average
      av = ((double)k * av + d) / (double)(k + 1);

      // get values of response variable
      yl[k] = data.at(smpl, resp);
    }

    // number of samples in current node
    nodecnt = rowMaskVec.size();

    // sum of all indep. variable values in current node
    sumnode = av * nodecnt;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }


      // create class indexes of current row
      idxRowVec.clear();
      for (j = 0; j < rowMaskVec.size(); ++j) {
        valCell = data.at(rowMaskVec[j], *idxCol);
        idxCell = data.find(valCell, data.varClassIndexer[*idxCol]);
        idxRowVec.push_back(idxCell);
      }

      // get variable's classes
      Helper::getClasses<uli_t>(idxRowVec, classVarVec);

      // break if there is just one class
      if (classVarVec.size() == 1) {
        ++non;

        if (non > 3 * data.par.ncol) { // too many variables were pure, so stop it
          fitFctPar.stopGrowing = true;
          return classifier;
        }

        ++idxCol;
        continue;
      }

      // check all combinations
      // => check all sets in power set 2^(numberofclasses-1) - 1
      // (avoiding duplicates)

      // get powerset (binary representation)
      powerSet.reset();
      powerSet.resize(classVarVec.size(), false);

      // check all sets
      do {

	      // init sums and number of samples
	      suml = npopl = 0;
	      sumr = sumnode;
	      npopr = nodecnt;

				// get the differences
				idxRow = rowMaskVec.begin();
				j = 0;
				while (idxRow != rowMaskVec.end()) {

					if (powerSet[idxRowVec[j]]) {

						val = data.at(*idxRow, depVar);

						// adjust sums and sample counts
						suml += val;
						sumr -= val;
						++npopl;
						--npopr;

					}
					++j;
					++idxRow;
				}

				critvar = (pow(suml, 2) / npopl) + (pow(sumr, 2) / npopr);

				// current indep. variable had got the best split
				if (critvar > critmax) {
					Helper::getPowerSet(powerSet, classVarVec, setCur);
					critmax = critvar;
					sosBestIdx = *idxCol;
					bestSubSet = setCur;
				}

				// get next set
        Helper::inc(powerSet);

      } while (!powerSet[powerSet.size() - 1]);

      ++idxCol;
    }

    // output dec. of sos
    if (critmax == 0) {
      fitFctPar.stopGrowing = true;
      decOfSos = 0;
    } else {
			// store final decrease of sum of squares
			decOfSos = critmax - (pow(sumnode, 2) / nodecnt);

      // get output vector
      Helper::getIndexedVector(
          data.varClassIndexer[sosBestIdx], bestSubSet, colVec);

      // build classifier
      classifier->setVarID(sosBestIdx);
      classifier->setClassSet(colVec);
      classifier->setMissingCode(data.getMissingCode());

      //add importance
      ((TImportance<T> *) fitFctPar.importance)->add(
          classifier->getVarID(), decOfSos / data.par.ntree);

#ifdef FITTINGFCTDEBUG2
      *fitFctPar.io->outVerbose << "twoingBestIdx " << sosBestIdx << std::endl;
      *fitFctPar.io->outVerbose << "bestSubSet(indexes) ";
      Helper::printVec(bestSubSet);
      *fitFctPar.io->outVerbose << "decOfTwoing " << decOfSos << std::endl << std::endl;
      *fitFctPar.io->outVerbose << "bestSubSet ";
      Helper::printVec(colVec);
      *fitFctPar.io->outVerbose << "classifier " << *classifier << std::endl;
#endif
    }
    return classifier;
  }

  /*
   * CARTtwoing
   *
   * Find the best nominal variable with categorical outcome.
   * The twoing method is described in CART (Leo Breiman, et al. 1983).
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTtwoing(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    //OUTPUT
    double &decOfTwoing = *fitFctPar.treeImprovement;
    uli_t &twoingBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    uli_t idxCell;
    T valCell;
    std::vector<T> colVec;
    std::vector<uli_t> setCur;
    std::vector<uli_t> bestSubSet;
    std::vector<uli_t> classVarVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> tmpclass;
    std::vector<uli_t> tclasspop;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxRowVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    SClassAtom<T> *classifier = new SClassAtom<T> ();
    boost::dynamic_bitset<> powerSet;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    double twoing;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    //double stdDecOfTwoing = 0.0;
    double ntwmofbest;
    double pln, pld, prn, prd, pno;

    a = b = ad = bd = 0.0;
    decOfTwoing = -std::numeric_limits<double>::infinity();

    // initial count

    k = j = 0;
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit = nl;
    idxVec.clear();
    idxRow = rowMaskVec.begin();
    j = 0;
    Nt = pno = 0;
    while (idxRow != rowMaskVec.end()) {
      // get index for next class in y
      idxVec.push_back(data.depIdxVec[*idxRow]);

      // number of samples (weighted) of node
      nrinit[idxVec[j]] += rowWeight[*idxRow];
      Nt += rowWeight[*idxRow];

      ++j;
      ++idxRow;
    }
    ntwmofbest = Nt;

    for (j = 0; j < nrinit.size(); ++j)
      pno += nrinit[j] * nrinit[j];

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      twoing = 0.0;

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "var " << data.varNames[*idxCol] << std::endl;
      *fitFctPar.io->outVerbose << "column ";
#endif
      // create class indexes of current row
      idxRowVec.clear();
      for (j = 0; j < rowMaskVec.size(); ++j) {
        valCell = data.at(rowMaskVec[j], *idxCol);
#ifdef FITTINGFCTDEBUG
        *fitFctPar.io->outVerbose << valCell << " ";
#endif
        idxCell = data.find(valCell, data.varClassIndexer[*idxCol]);
        idxRowVec.push_back(idxCell);
      }
#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << std::endl;
#endif

      // get variable's classes
      Helper::getClasses<uli_t>(idxRowVec, classVarVec);

      // break if there is just one class
      if (classVarVec.size() == 1) {
        //twoing = stdDecOfTwoing;
        twoing = 0;
        if (twoing > decOfTwoing) {
          setCur.clear();
          decOfTwoing = twoing;
          twoingBestIdx = *idxCol;
          bestSubSet = setCur;
        }
        ++idxCol;
        continue;
      }

      // check all combinations
      // => check all sets in power set 2^(numberofclasses-1) - 1
      // (avoiding duplicates)

      // get powerset (binary representation)
      powerSet.reset();
      powerSet.resize(classVarVec.size(), false);

#ifdef FITTINGFCTDEBUG2
      *fitFctPar.io->outVerbose << "idxVec ";
      Helper::printVec<uli_t>(idxVec);
      *fitFctPar.io->outVerbose << "idxRowVec ";
      Helper::printVec<uli_t>(idxRowVec);
      *fitFctPar.io->outVerbose << "classVarVec ";
      Helper::printVec<uli_t>(classVarVec);
#endif

      // check all sets
      do {
        // count number of samples per class in y in sub set
#ifdef __DEBUG__
			*fitFctPar.io->outVerbose << "tmpclass.resize(data.varClassIndexer[data.par.depVar].size(), 0) = " <<
				"tmpclass.resize(data.varClassIndexer[" << data.par.depVar << "].size(), 0)" <<
				"tmpclass.resize(" << data.varClassIndexer[data.par.depVar].size() << ", 0)" << std::endl;
			*fitFctPar.io->outVerbose << "data.varClassIndexer = " << data.varClassIndexer << std::endl;
#endif	
        tmpclass.clear();
        tmpclass.resize(data.varClassIndexer[data.par.depVar].size(), 0);
        for (j = 0; j < idxVec.size(); ++j)
          if (powerSet[idxRowVec[j]])
            tmpclass[idxVec[j]] += (uli_t) rowWeight[rowMaskVec[j]];

#ifdef FITTINGFCTDEBUG2
        *fitFctPar.io->outVerbose << "powerSet " << powerSet << std::endl;
        *fitFctPar.io->outVerbose << "tmpclass ";
        Helper::printVec(tmpclass);
        *fitFctPar.io->outVerbose << "nrinit ";
        Helper::printVec(nrinit);
#endif

        prn = pln = pld = prd = 0.0;

        for (j = 0; j < tmpclass.size(); ++j) {
          pln += tmpclass[j] * tmpclass[j];
          pld += tmpclass[j];
        }

        for (j = 0; j < tmpclass.size(); ++j) {
          tmpclass[j] = (uli_t) nrinit[j] - tmpclass[j];
          prn += tmpclass[j] * tmpclass[j];
        }

        prd = (Nt - pld);

        a = (pld == 0) ? (0)
            : (pln / pld);
        b = (prd == 0) ? (0)
            : (prn / prd);

        twoing = a + b;

#ifdef FITTINGFCTDEBUG
        *fitFctPar.io->outVerbose << "Nt " << Nt << std::endl;
        *fitFctPar.io->outVerbose << "pld " << pld << std::endl;
        *fitFctPar.io->outVerbose << "prn " << prn << std::endl;
        *fitFctPar.io->outVerbose << "pln " << pln << std::endl;
        *fitFctPar.io->outVerbose << "pno " << pno << std::endl;
        *fitFctPar.io->outVerbose << "twoing " << twoing << std::endl << std::endl;
#endif

        // is it a good twoing?
        if (twoing > decOfTwoing) {
          Helper::getPowerSet(powerSet, classVarVec, setCur);
          decOfTwoing = twoing;
          twoingBestIdx = *idxCol;
          bestSubSet = setCur;
        }

        // get next set
        Helper::inc(powerSet);

      } while (!powerSet[powerSet.size() - 1]);

      ++idxCol;
    }

    decOfTwoing = decOfTwoing - pno / Nt;

    // get output vector
    Helper::getIndexedVector(
        data.varClassIndexer[twoingBestIdx], bestSubSet, colVec);

    // build classifier
    classifier->setVarID(twoingBestIdx);
    classifier->setClassSet(colVec);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T>*) fitFctPar.importance)->add(twoingBestIdx, decOfTwoing
        / data.par.ntree);

#ifdef FITTINGFCTDEBUG2
    *fitFctPar.io->outVerbose << "twoingBestIdx " << twoingBestIdx << std::endl;
    *fitFctPar.io->outVerbose << "bestSubSet(indexes) ";
    Helper::printVec(bestSubSet);
    *fitFctPar.io->outVerbose << "decOfTwoing " << decOfTwoing << std::endl << std::endl;
    *fitFctPar.io->outVerbose << "bestSubSet ";
    Helper::printVec(colVec);
    *fitFctPar.io->outVerbose << "classifier " << *classifier << std::endl;
#endif

    return classifier;
  }

  /*
   * CART two way interaction, additive model
   *
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTtwoWayIntAdd(FittingFctPar<T> &fitFctPar) {
    IAMClassAtom<T> *classifier = new IAMClassAtom<T> ();

    // get all paired combinations of selected variables

    // process each pair

    // merge both columns of the variable pair:
    // merge method: a + b (Addition)

    // do a typical gini search on merged column

    // save best variable pair

    // create classifier object


    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    //OUTPUT
    double &decOfTwoing = *fitFctPar.treeImprovement;
    uli_t &twoingBestIdx1 = *fitFctPar.bestVar;
    uli_t twoingBestIdx2;

    uli_t j, k;
    uli_t idxCell1, idxCell2, idxCount;
    T valCell;
    std::vector<T> colVec;
    std::vector<uli_t> setCur;
    std::vector<uli_t> bestSubSet;
    std::vector<uli_t> classVarVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> tmpclass;
    std::vector<uli_t> tclasspop;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxRowVec;
    std::vector<std::vector<uli_t> > idxMapMat;
    std::vector<std::pair<uli_t, uli_t> > idxBackMapMat;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    uli_t var1, var2;
    boost::dynamic_bitset<> powerSet;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    double twoing;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    //double stdDecOfTwoing = 0.0;
    double ntwmofbest;
    double pln, pld, prn, prd, pno;

    var1 = var2 = idxCell1 = idxCell2 = 0;

    a = b = ad = bd = 0.0;
    decOfTwoing = -std::numeric_limits<double>::infinity();

    // initial count

    k = j = 0;
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);
    idxVec.clear();
    idxRow = rowMaskVec.begin();
    j = 0;
    Nt = pno = 0;
    while (idxRow != rowMaskVec.end()) {
      // get index for next class in y
      idxVec.push_back(data.depIdxVec[*idxRow]);

      // number of samples (weighted) of node
      nrinit[idxVec[j]] += rowWeight[*idxRow];
      Nt += rowWeight[*idxRow];

      ++j;
      ++idxRow;
    }
    ntwmofbest = Nt;

    for (j = 0; j < nrinit.size(); ++j)
      pno += nrinit[j] * nrinit[j];

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      var1 = *idxCol;
      if ((idxCol + 1) == colMaskVec.end())
        break;
      var2 = *(idxCol + 1);

      twoing = 0.0;

      // create class indexes of current row
      idxCount = 0;
      idxRowVec.clear();
      idxMapMat = std::vector<std::vector<uli_t> >(
          data.varClassIndexer[var1].size(), std::vector<uli_t>(
              data.varClassIndexer[var2].size(),
              std::numeric_limits<uli_t>::max()));

      for (j = 0; j < rowMaskVec.size(); ++j) {
        // get indexes of variable 1 and 2
        valCell = data.at(rowMaskVec[j], var1);
        idxCell1 = data.find(valCell, data.varClassIndexer[var1]);

        valCell = data.at(rowMaskVec[j], var2);
        idxCell2 = data.find(valCell, data.varClassIndexer[var2]);

        // map two indexes to one new index and store the mapping
        if (idxMapMat[idxCell1][idxCell2] == std::numeric_limits<uli_t>::max()) {
          idxMapMat[idxCell1][idxCell2] = idxCount;
          idxBackMapMat.push_back(std::make_pair(idxCell1, idxCell2));
          ++idxCount;
        }

        idxRowVec.push_back(idxMapMat[idxCell1][idxCell2]);
      }

      // get variable's classes
      Helper::getClasses<uli_t>(idxRowVec, classVarVec);

      // break if there is just one class
      if (classVarVec.size() == 1) {
        //twoing = stdDecOfTwoing;
        twoing = 0;
        if (twoing > decOfTwoing) {
          setCur.clear();
          decOfTwoing = twoing;
          twoingBestIdx1 = var1;
          twoingBestIdx2 = var2;
          bestSubSet = setCur;
        }
        ++idxCol;
        continue;
      }

      // check all combinations
      // => check all sets in power set 2^(numberofclasses-1) - 1
      // (avoiding duplicates)

      // get powerset (binary representation)
      powerSet.reset();
      powerSet.resize(classVarVec.size(), false);

      // check all sets
      do {
        // count number of samples per class in y in sub set
        tmpclass.clear();
        tmpclass.resize(data.varClassIndexer[data.par.depVar].size(), 0);
        for (j = 0; j < idxVec.size(); ++j)
          if (powerSet[idxRowVec[j]])
            tmpclass[idxVec[j]] += (uli_t) rowWeight[rowMaskVec[j]];

        prn = pln = pld = prd = 0.0;

        for (j = 0; j < tmpclass.size(); ++j) {
          pln += tmpclass[j] * tmpclass[j];
          pld += tmpclass[j];
        }

        for (j = 0; j < tmpclass.size(); ++j) {
          tmpclass[j] = (uli_t) nrinit[j] - tmpclass[j];
          prn += tmpclass[j] * tmpclass[j];
        }

        prd = (Nt - pld);

        a = (pld == 0) ? (0)
            : (pln / pld);
        b = (prd == 0) ? (0)
            : (prn / prd);

        twoing = a + b;

        // is it a good twoing?
        if (twoing > decOfTwoing) {
          Helper::getPowerSet(powerSet, classVarVec, setCur);
          decOfTwoing = twoing;
          twoingBestIdx1 = var1;
          twoingBestIdx2 = var2;
          bestSubSet = setCur;
        }

        // get next set
        Helper::inc(powerSet);

      } while (!powerSet[powerSet.size() - 1]);

      ++idxCol;
    }

    decOfTwoing -= pno / Nt;

    // get output vector
    //Helper::getIndexedVector(
    //  data.varClassIndexer[twoingBestIdx], bestSubSet, colVec);

    // build classifier
    classifier->varID.push_back(twoingBestIdx1);
    classifier->varID.push_back(twoingBestIdx2);
    classifier->setMissingCode(data.getMissingCode());

    // restore values of variable 1 and 2
    uli_t idx;
    classifier->classSet.push_back(std::vector<T>());
    classifier->classSet.push_back(std::vector<T>());
    for (uli_t i = 0; i < bestSubSet.size(); ++i) {
      idx = idxBackMapMat[bestSubSet[i]].first;
      classifier->classSet[0].push_back(
          data.varClassIndexer[twoingBestIdx1][idx]);

      idx = idxBackMapMat[bestSubSet[i]].second;
      classifier->classSet[1].push_back(
          data.varClassIndexer[twoingBestIdx2][idx]);
    }

    //add importance
    ((IAM2WayImportance<T>*) fitFctPar.importance)->add(
        twoingBestIdx1, twoingBestIdx2, decOfTwoing / data.par.ntree);

    return classifier;
  }


  /*
   * check for sse
   */
	template<class T>
		static void checkUsabilitySse(DataFrame<T> &data) {
			bool hasSSE;

			// we have to use local numbered labels in asm code
			// because optimize might duplicate code
			// 0f: jump forward to label 0
			// 0b: jump backward to label 0
			__asm__ __volatile__(
					"movl $1, %%eax\n\t"
					"cpuid\n\t"
					"test $0x002000000, %%edx\n\t"
					"jnz 0f\n\t"
					"movl $0, %%eax\n\t"
					"jmp 1f\n\t"
					"0:\n\t"
					"movl $1, %%eax\n\t"
					"1:"
					:"=a"(hasSSE)
					:   
					:); 
			if (!hasSSE) {
				throw Exception(ERRORCODE_63);
			}
		}

	
  /*
   * CARTgini
   *
   * Find the best numerical/ordinal variable with categorical outcome.
   * The gini method is described in CART (Leo Breiman, et al. 1984).
   * It regards missing values.
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTgini(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<double> nl;
    std::vector<double> nr;
    std::vector<double> nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff = 0;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    T cutoff;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini = 0.0;
    double ntwm, ntwmofbest;

    a = b = ad = bd = 0.0;
    decOfGini = -std::numeric_limits<double>::infinity();

    // initial count

    k = j = 0;
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);
    idxVec.clear();
    idxRow = rowMaskVec.begin();
    j = 0;
    Nt = 0;
    while (idxRow != rowMaskVec.end()) {
      idxVec.push_back(data.depIdxVec[*idxRow]);
      nrinit[idxVec[j]] += rowWeight[*idxRow];
      ++j;
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    b = 0.0;
    bd = 0.0;
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;
    Nt = bd;
    ntwmofbest = Nt;
    ntwm = Nt;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      gini = 0.0;

      // get variable's classes
      data.getRowMaskedCol(*idxCol, rowMaskVec, colVec);
      Helper::getClasses<T>(colVec, classVarVec);
      sort(classVarVec.begin(), classVarVec.end());
      //data.getRowMaskedClassesInCol(*idxCol, rowMaskVec, classVarVec);

      // break if there is just one class
      if (classVarVec.size() == 1) {
        gini = stdDecOfGini;
        if (gini > decOfGini) {
          decOfGini = gini;
          giniBestIdx = *idxCol;
          bestCutoff = (double) classVarVec[0];
        }
        ++idxCol;
        continue;
      }

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }
      nr = nrinit;

      // iterate over all cutoffs
      for (j = 0; j < classVarVec.size() - 1; ++j) {
        // intrinsic cutoff between two value which is
        // (v_i + v_{i+1})/2
        cutoff = (classVarVec[j]);

        k = 0;
        idxRow = rowMaskVec.begin();
        // get the differences
        while (idxRow != rowMaskVec.end()) {
          if (data.at(*idxRow, *idxCol) == cutoff) {
            nl[idxVec[k]] += rowWeight[*idxRow];
            nr[idxVec[k]] -= rowWeight[*idxRow];
          }
          ++k;
          ++idxRow;
        }

        // calculate updated gini
        // left child.
        itlr = nl.begin();
        a = 0.0;
        ad = 0.0;
        while (itlr != nl.end()) {
          a += (*itlr) * (*itlr);
          ad += (*itlr);
          ++itlr;
        }
        (ad != 0.0) ? (a /= ad)
            : (a = 0.0);

        // right child.
        itlr = nr.begin();
        b = 0.0;
        bd = 0.0;
        while (itlr != nr.end()) {
          b += (*itlr) * (*itlr);
          bd += (*itlr);
          ++itlr;
        }
        (bd != 0.0) ? (b /= bd)
            : (b = 0.0);

        gini = a + b;

#ifdef FITTINGFCTDEBUG2
        Helper::printVec<double>(nl);
        Helper::printVec<double>(nr);
        Helper::printVec<double>(rowWeight);
        Helper::printVec<T>(colVec);
        Helper::printVec<T>(classVarVec);
        Helper::printVec<uli_t>(rowMaskVec);
        *fitFctPar.io->outVerbose << (double) data.at(0, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << (double) data.at(3, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << (double) data.at(4, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << (double) data.at(8, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
        *fitFctPar.io->outVerbose << "gini " << gini << std::endl;
        *fitFctPar.io->outVerbose << "cutoff " << cutoff << std::endl;
        *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
        *fitFctPar.io->outVerbose << "giniBestIdx " << giniBestIdx << std::endl;
        *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
        *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif
        // is it a good gini?
        if ((gini > decOfGini) ||
						((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
          decOfGini = gini;
          giniBestIdx = *idxCol;
          bestCutoff = ((double) cutoff + (double) classVarVec[j + 1]) / 2.0;
          ntwmofbest = ntwm;
        }

      }
      ++idxCol;
    }
#ifdef FITTINGFCTDEBUG
    *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
    *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
    *fitFctPar.io->outVerbose << "nodeGini " << nodeGini << std::endl;
    *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
    *fitFctPar.io->outVerbose << "(T)bestCutoff " << (T) bestCutoff << std::endl;
#endif

    // weighting with sample size
    // what rf5 of breiman/cutler do
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // calc. the true dec. of gini.
    //decOfGini = nodeGini - 1 + decOfGini / ntwmofbest;

#ifdef FITTINGFCTDEBUG
    *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif
    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T>*) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

  template<class T>
  static ClassAtom<T, uli_t> *CARTginiNoDenominator(FittingFctPar<T> &fitFctPar) {
    INode<T> *&parent = fitFctPar.parent;
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    uli_t iteration = fitFctPar.iteration;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> leftOutcomeVec, rightOutcomeVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff = 0;
    double a, b, ad, bd, bestad, bestbd, Nt;
    double cutoff;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini = 0.0;
    double ntwm, ntwmofbest;

    a = b = ad = bd = bestad = bestbd = 0.0;
    decOfGini = -std::numeric_limits<double>::infinity();

    // initial count

    k = j = 0;
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);
    idxVec.clear();
    idxRow = rowMaskVec.begin();
    j = 0;
    Nt = 0;
    while (idxRow != rowMaskVec.end()) {
      idxVec.push_back(data.depIdxVec[*idxRow]);
      nrinit[idxVec[j]] += rowWeight[*idxRow];
      Nt += rowWeight[*idxRow];
      ++j;
      ++idxRow;
    }
    ntwmofbest = Nt;

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    b = 0.0;
    bd = 0.0;
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;
    bestbd = bd;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      gini = 0.0;

      // get variable's classes
      data.getRowMaskedCol(*idxCol, rowMaskVec, colVec);
      Helper::getClasses<T>(colVec, classVarVec);
      sort(classVarVec.begin(), classVarVec.end());
      //data.getRowMaskedClassesInCol(*idxCol, rowMaskVec, classVarVec);

      // break if there is just one class
      if (classVarVec.size() == 1) {
        gini = b;
        if (gini > decOfGini) {
          decOfGini = gini;
          giniBestIdx = *idxCol;
          bestCutoff = (double) classVarVec[0];
        }
        ++idxCol;
        continue;
      }

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }
      nr = nrinit;
      ntwm = Nt; //"degrade" missing samples

      // iterate over all cutoffs
      for (j = 0; j < classVarVec.size() - 1; ++j) {
        // intrinsic cutoff between two value which is
        // (v_i + v_{i+1})/2
        cutoff = (double) (classVarVec[j]);

        k = 0;
        idxRow = rowMaskVec.begin();
        // get the differences
        while (idxRow != rowMaskVec.end()) {
          if (data.at(*idxRow, *idxCol) == cutoff) {
            nl[idxVec[k]] += rowWeight[*idxRow];
            nr[idxVec[k]] -= rowWeight[*idxRow];
          }
          ++k;
          ++idxRow;
        }

        // calculate updated gini
        // left child.
        itlr = nl.begin();
        a = 0.0;
        ad = 0.0;
        while (itlr != nl.end()) {
          a += (*itlr) * (*itlr);
          ad += (*itlr);
          ++itlr;
        }

        // right child.
        itlr = nr.begin();
        b = 0.0;
        bd = 0.0;
        while (itlr != nr.end()) {
          b += (*itlr) * (*itlr);
          bd += (*itlr);
          ++itlr;
        }

        (ad != 0.0) ? (b *= ad)
            : (a = 0.0);
        (bd != 0.0) ? (a *= bd)
            : (b = 0.0);

        gini = a + b;

#ifdef FITTINGFCTDEBUG2
        Helper::printVec<double>(nl);
        Helper::printVec<double>(nr);
        Helper::printVec<double>(rowWeight);
        Helper::printVec<T>(colVec);
        Helper::printVec<T>(classVarVec);
        Helper::printVec<uli_t>(rowMaskVec);
        *fitFctPar.io->outVerbose << data.at(0, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << data.at(3, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << data.at(4, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << data.at(8, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
        *fitFctPar.io->outVerbose << "gini " << gini << std::endl;
        *fitFctPar.io->outVerbose << "cutoff " << cutoff << std::endl;
        *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
        *fitFctPar.io->outVerbose << "giniBestIdx " << giniBestIdx << std::endl;
        *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
        *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif

        // is it a good gini?
        if ((gini > decOfGini) ||
						((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
          decOfGini = gini;
          giniBestIdx = *idxCol;
          bestCutoff = (cutoff + (double) classVarVec[j + 1]) / 2.0;
          bestad = ad;
          bestbd = bd;
          ntwmofbest = ntwm;
        }

      }
      ++idxCol;
    }

    *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
    *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
    *fitFctPar.io->outVerbose << "nodeGini " << nodeGini << std::endl;
    *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
    *fitFctPar.io->outVerbose << "(T)bestCutoff " << (T) bestCutoff << std::endl;

    // output dec. of gini
    decOfGini /= ntwmofbest;
    if (bestad != 0.0)
      decOfGini /= bestad;
    if (bestbd != 0.0)
      decOfGini /= bestbd;
    decOfGini += nodeGini - 1;

    *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;

    // build classifier
    classifier->setVarID(giniBestIdx);
    *fitFctPar.io->outVerbose << "1" << std::endl;
    classifier->setThreshold((T) bestCutoff);
    *fitFctPar.io->outVerbose << "2" << std::endl;
    classifier->setMissingCode(data.getMissingCode());
    *fitFctPar.io->outVerbose << "3" << std::endl;

    //add importance
    fitFctPar.importance->add(classifier->getVarID(), decOfGini);
    *fitFctPar.io->outVerbose << "4" << std::endl;

    return classifier;
  }

  /*
   * CARTgini
   *
   * Find the best numerical/ordinal variable with categorical outcome.
   * The gini method is described in CART (Leo Breiman, et al. 1983).
   * It regards missing values.
   */
  /*
   template <class T >
   static ClassAtom<T, uli_t > *CARTginiWithMissings(//INPUT
   INode<T > *parent,
   DataFrame<T > &data,
   uli_t depVar,
   std::vector<uli_t > &rowMaskVec,
   std::vector<uli_t > &colMaskVec,
   std::vector<double > &rowWeight,
   double nodeGini,
   uli_t iteration,
   //OUTPUT
   double &decOfGini,
   uli_t &giniBestIdx
   ) {

   uli_t j, k;
   std::vector<T > colVec;
   std::vector<T > classVarVec;
   std::vector<T > giniVec;
   std::vector<T > outcomeVec;
   std::vector<double > nl, nr, nrinit;
   std::vector<double >::iterator itlr;
   std::vector<uli_t > idxVec;
   std::vector<uli_t > maskVec;
   std::vector<uli_t > classMaskVec;
   std::vector<uli_t > rowMaskForDepVec;
   std::vector<uli_t > classOutcomeVec;
   typename std::vector<uli_t >::iterator idxRow;
   typename std::vector<uli_t >::iterator idxCol;
   T2ClassAtom<T > *classifier = new T2ClassAtom<T >();
   double bestCutoff = 0;
   double a, b, ad, bd, Nt; // Nt = number of samples in current node
   double cutoff;
   double gini;
   uli_t depVarClassSize = data.varCategories[depVar].size();
   double stdDecOfGini = 0.0;
   double ntwm, ntwmofbest;

   a = b = ad = bd = 0.0;
   decOfGini = -std::numeric_limits<double>::infinity();
   T missingcode = data.getMissingCode();

   // initial count

   k = j = 0;
   nl.clear();
   nl.resize(depVarClassSize, 0.0);
   nrinit.clear();
   nrinit.resize(depVarClassSize, 0.0);
   idxVec.clear();
   idxRow = rowMaskVec.begin();
   j = 0;
   Nt = 0;
   while (idxRow != rowMaskVec.end()) {
   idxVec.push_back(data.depIdxVec[*idxRow]);
   nrinit[idxVec[j]] += rowWeight[*idxRow];
   Nt += rowWeight[*idxRow];
   ++j;
   ++idxRow;
   }
   ntwmofbest = Nt;

   // 1 + gini (without normalisation) of the whole node
   itlr = nrinit.begin();
   b = 0.0;
   bd = 0.0;
   while (itlr != nrinit.end()) {
   b  += (*itlr) * (*itlr);
   bd += (*itlr);
   ++itlr;
   }
   (bd != 0.0)?(b /= bd):(b = 0.0);
   stdDecOfGini = b;

   idxCol = colMaskVec.begin();
   while (idxCol != colMaskVec.end()) {
   //skip dep var
   if (*idxCol == depVar) {
   ++idxCol;
   continue;
   }

   gini = 0.0;

   // get variable's classes
   data.getRowMaskedCol(*idxCol, rowMaskVec, colVec);
   Helper::getClasses<T >(colVec, classVarVec);
   sort(classVarVec.begin(), classVarVec.end());

   // break if there is just one class
   if (classVarVec.size() == 1 && classVarVec[0] == missingcode) {
   ++idxCol;
   continue;
   } else if (
   (classVarVec.size() <= 2 && classVarVec[0] == missingcode) ||
   (classVarVec.size() <= 2 && classVarVec[1] == missingcode) ||
   classVarVec.size() == 1) {
   gini = stdDecOfGini;
   if (gini > decOfGini) {
   decOfGini = gini;
   giniBestIdx = *idxCol;
   if (classVarVec.size() == 1 ||
   (classVarVec.size() <= 2 && classVarVec[1] == missingcode))
   bestCutoff = (double)classVarVec[0];
   else
   bestCutoff = (double)classVarVec[1];
   }
   ++idxCol;
   continue;
   }

   // initial reset of the left node
   itlr = nl.begin();
   while (itlr != nl.end()) {
   *itlr = 0.0;
   ++itlr;
   }
   nr = nrinit;
   ntwm = Nt; //"degrade" missing samples

   // iterate over all cutoffs
   for (j = 0; j < classVarVec.size() - 1; ++j) {
   // intrinsic cutoff between two value which is
   // (v_i + v_{i+1})/2
   cutoff = (double)(classVarVec[j]);

   k = 0;
   idxRow = rowMaskVec.begin();
   if (classVarVec[j] == missingcode) {
   // degrade missings
   while (idxRow != rowMaskVec.end()) {
   if (data.at(*idxRow, *idxCol) == missingcode) {
   ntwm -= rowWeight[*idxRow];
   }
   ++k;
   ++idxRow;
   }
   continue;
   } else {
   // get the differences
   while (idxRow != rowMaskVec.end()) {
   if (data.at(*idxRow, *idxCol) == cutoff) {
   nl[idxVec[k]] += rowWeight[*idxRow];
   nr[idxVec[k]] -= rowWeight[*idxRow];
   }
   ++k;
   ++idxRow;
   }
   }

   // calculate updated gini
   // left child.
   itlr = nl.begin();
   a = 0.0;
   ad = 0.0;
   while (itlr != nl.end()) {
   a  += (*itlr) * (*itlr);
   ad += (*itlr);
   ++itlr;
   }
   (ad != 0.0)?(a /= ad):(a = 0.0);

   // right child.
   itlr = nr.begin();
   b = 0.0;
   bd = 0.0;
   while (itlr != nr.end()) {
   b  += (*itlr) * (*itlr);
   bd += (*itlr);
   ++itlr;
   }
   (bd != 0.0)?(b /= bd):(b = 0.0);

   gini = a + b;

   // is it a good gini?
   if (gini > decOfGini) {
   decOfGini = gini;
   giniBestIdx = *idxCol;
   bestCutoff = (cutoff + (double)classVarVec[j+1]) / 2.0;
   ntwmofbest = ntwm;
   }
   }
   ++idxCol;
   }

   // output dec. of gini
   decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

   // build classifier
   classifier->setVarID(giniBestIdx);
   classifier->setThreshold((T)bestCutoff);
   classifier->setMissingCode(data.getMissingCode());

   //add importance
   fitFctPar.importance->add(classifier->getVarID(), decOfGini);

   return classifier;
   }
   */

  /*
   * CARTginiSrt
   *
   * Find the best numerical/ordinal variable with categorical outcome.
   * The gini method is described in CART (Leo Breiman, et al. 1983).
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTginiSrtOld(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxMissVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    T cutoffT, cutoffTlast;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini;
    double ntwm, ntwmofbest;

    a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;
    k = j = 0;
    cutoffT = cutoffTlast = 0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);

    // init count
    idxRow = rowMaskVec.begin();
    while (idxRow != rowMaskVec.end()) {
      nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    Nt = bd;
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "*idxCol.. " << *idxCol << std::endl;
#endif
      gini = 0.0;

      // get column data information
			data.getRowMaskedColNoMissing2AndPair(*idxCol, rowMaskVec, rowPairVec, idxMissVec);

#ifdef FITTINGFCTDEBUG2
      *fitFctPar.io->outVerbose << "rowPairVec.size() " << rowPairVec.size() << std::endl;
      Helper::printVec(colVec);
#endif

      // if data contains only missings then skip the current column
      if (rowPairVec.size() == 0) {
        ++idxCol;
        continue;
      }

      // decount missings
      nr = nrinit;
      ntwmofbest = ntwm = Nt;

      sort(rowPairVec.begin(), rowPairVec.end());

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }

      // iterate over all cutoffs in rowPairVec
      cutoffT = missingcode;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        // change cutoff if necessary
        // and calc gini index
        if (cutoffT != itRowPairVec->first) {
          // first iteration check
          (missingcode == cutoffT) ? cutoffTlast = itRowPairVec->first
              : cutoffTlast = cutoffT;
          cutoffT = itRowPairVec->first;

          // calculate updated gini
          // left child.
          itlr = nl.begin();
          a = 0.0;
          ad = 0.0;
          while (itlr != nl.end()) {
            a += (*itlr) * (*itlr);
            ad += (*itlr);
            ++itlr;
          }
          (ad != 0.0) ? (a /= ad)
              : (a = 0.0);

          // right child.
          itlr = nr.begin();
          b = 0.0;
          bd = 0.0;
          while (itlr != nr.end()) {
            b += (*itlr) * (*itlr);
            bd += (*itlr);
            ++itlr;
          }
          (bd != 0.0) ? (b /= bd)
              : (b = 0.0);

          gini = a + b;

#ifdef FITTINGFCTDEBUG2
          Helper::printVec<double>(nl);
          Helper::printVec<double>(nr);
          Helper::printVec<double>(rowWeight);
          Helper::printVec<T>(colVec);
          Helper::printVec<T>(classVarVec);
          Helper::printVec<uli_t>(rowMaskVec);
          *fitFctPar.io->outVerbose << (double) data.at(0, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(3, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(4, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(8, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
          *fitFctPar.io->outVerbose << "gini " << gini << std::endl;
          *fitFctPar.io->outVerbose << "cutoffT " << cutoffT << std::endl;
          *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
          *fitFctPar.io->outVerbose << "giniBestIdx " << giniBestIdx << std::endl;
          *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
          *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif

          // is it a good gini check
					if ((gini > decOfGini) ||
							((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
						decOfGini = gini;
						giniBestIdx = *idxCol;
						bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
						ntwmofbest = ntwm;
					}
        }

#ifdef FITTINGFCTDEBUG
        *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
        *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
        *fitFctPar.io->outVerbose << "nodeGini " << nodeGini << std::endl;
        *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
        *fitFctPar.io->outVerbose << "(T)bestCutoff " << (T) bestCutoff << std::endl;
#endif

        // get the differences
        k = data.depIdxVec[itRowPairVec->second];
        nl[k] += rowWeight[itRowPairVec->second];
        nr[k] -= rowWeight[itRowPairVec->second];

        ++itRowPairVec;
      }
      ++idxCol;
    }

    // weighting with sample size
    // what rf5 of breiman/cutler do
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // calc. the true dec. of gini.
    //decOfGini = nodeGini - 1 + decOfGini / ntwmofbest;


    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

  /*
   * CARTginiSrt
   *
   * Find the best numerical/ordinal variable with categorical outcome.
   * The gini method is described in CART (Leo Breiman, et al. 1983).
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTginiSrt(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxMissVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff;
    double w, a, b, ad, bd, Nt, Nt2, noml, denoml, nomr, denomr; // Nt = number of samples in current node
    T cutoffT, cutoffTlast;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini;
    double ntwm, ntwmofbest;

    w = a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;

    k = j = 0;
    cutoffT = cutoffTlast = 0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);

    // init count
    idxRow = rowMaskVec.begin();
    while (idxRow != rowMaskVec.end()) {
      nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    Nt = bd;
		Nt2 = b;
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "*idxCol.. " << *idxCol << std::endl;
#endif
      gini = 0.0;

      // get column data information
			data.getRowMaskedColNoMissing2AndPair(*idxCol, rowMaskVec, rowPairVec, idxMissVec);

#ifdef FITTINGFCTDEBUG2
      *fitFctPar.io->outVerbose << "rowPairVec.size() " << rowPairVec.size() << std::endl;
      Helper::printVec(colVec);
#endif

      // if data contains only missings then skip the current column
      if (rowPairVec.size() == 0) {
        ++idxCol;
        continue;
      }

      // decount missings
      nr = nrinit;
      ntwmofbest = ntwm = Nt;

      sort(rowPairVec.begin(), rowPairVec.end());

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }

		  noml = 0.0;
  		denoml = 0.0;
		  nomr = Nt2;
  		denomr = Nt;

      // iterate over all cutoffs in rowPairVec
      cutoffT = missingcode;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        // change cutoff if necessary
        // and calc gini index
        if (cutoffT != itRowPairVec->first) {
          // first iteration check
          (missingcode == cutoffT) ? cutoffTlast = itRowPairVec->first
              : cutoffTlast = cutoffT;
          cutoffT = itRowPairVec->first;

          gini = (noml / denoml) + (nomr / denomr);

#ifdef FITTINGFCTDEBUG2
          Helper::printVec<double>(nl);
          Helper::printVec<double>(nr);
          Helper::printVec<double>(rowWeight);
          Helper::printVec<T>(colVec);
          Helper::printVec<T>(classVarVec);
          Helper::printVec<uli_t>(rowMaskVec);
          *fitFctPar.io->outVerbose << (double) data.at(0, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(3, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(4, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(8, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
          *fitFctPar.io->outVerbose << "gini " << gini << std::endl;
          *fitFctPar.io->outVerbose << "cutoffT " << cutoffT << std::endl;
          *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
          *fitFctPar.io->outVerbose << "giniBestIdx " << giniBestIdx << std::endl;
          *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
          *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif

          // is it a good gini check
					if ((gini > decOfGini) ||
							((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
						decOfGini = gini;
						giniBestIdx = *idxCol;
						bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
						ntwmofbest = ntwm;
					}
        }

#ifdef FITTINGFCTDEBUG
        *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
        *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
        *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
        *fitFctPar.io->outVerbose << "(T)bestCutoff " << (T) bestCutoff << std::endl;
#endif

        // get the differences
        k = data.depIdxVec[itRowPairVec->second];
				w = rowWeight[itRowPairVec->second];

				// update contibution using binomial theorem
			  noml += w * (2 * nl[k] + w);
			  nomr += w * (-2 * nr[k] + w);

				// update denominators
  			denoml += w;
  			denomr -= w;

				// update class weights
				nl[k] += w;
				nr[k] -= w;

        ++itRowPairVec;
      }
      ++idxCol;
    }

    // weighting with sample size
    // what rf5 of breiman/cutler do
		decOfGini = decOfGini - Nt2/Nt;

    // calc. the true dec. of gini.


    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

	 /** 
		* \brief Cart srt using sse
		* 
		* @param fitFctPar
		* 
		* @return 
		*/
  template<class T>
  static ClassAtom<T, uli_t> *CARTginisse(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxMissVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    T cutoffT, cutoffTlast;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini;
    double ntwm, ntwmofbest;

    a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;
    k = j = 0;
    cutoffT = cutoffTlast = 0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);

    // init count
    idxRow = rowMaskVec.begin();
    while (idxRow != rowMaskVec.end()) {
      nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    Nt = bd;
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "*idxCol.. " << *idxCol << std::endl;
#endif
      gini = 0.0;

      // get column data information
			data.getRowMaskedColNoMissing2AndPair(*idxCol, rowMaskVec, rowPairVec, idxMissVec);
      //data.getRowMaskedColNoMissing2(*idxCol, rowMaskVec, colVec, idxVec, idxMissVec);

#ifdef FITTINGFCTDEBUG2
      *fitFctPar.io->outVerbose << "rowPairVec.size() " << rowPairVec.size() << std::endl;
      Helper::printVec(colVec);
#endif

      // if data contains only missings then skip the current column
      if (rowPairVec.size() == 0) {
        ++idxCol;
        continue;
      }

      // decount missings
      nr = nrinit;
      ntwmofbest = ntwm = Nt;

			sort(rowPairVec.begin(), rowPairVec.end());

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }

      // iterate over all cutoffs in rowPairVec
      cutoffT = missingcode;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        // change cutoff if necessary
        // and calc gini index
        if (cutoffT != itRowPairVec->first) {
          // first iteration check
          cutoffTlast = (missingcode == cutoffT) ? itRowPairVec->first : cutoffT;
          cutoffT = itRowPairVec->first;

          // calculate updated gini
          // left child.
          itlr = nl.begin();
          a = 0.0;
          ad = 0.0;
          while (itlr != nl.end()) {
            a += (*itlr) * (*itlr);
            ad += (*itlr);
            ++itlr;
          }
          (ad != 0.0) ? (a /= ad)
              : (a = 0.0);

          // right child.
          itlr = nr.begin();
          b = 0.0;
          bd = 0.0;
          while (itlr != nr.end()) {
            b += (*itlr) * (*itlr);
            bd += (*itlr);
            ++itlr;
          }
          (bd != 0.0) ? (b /= bd)
              : (b = 0.0);

          gini = a + b;

          // is it a good gini check
					if ((gini > decOfGini) ||
							((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
						decOfGini = gini;
						giniBestIdx = *idxCol;
            bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
            ntwmofbest = ntwm;
          }
        }

        // get the differences
        k = data.depIdxVec[itRowPairVec->second];
        nl[k] += rowWeight[itRowPairVec->second];
        nr[k] -= rowWeight[itRowPairVec->second];

        ++itRowPairVec;
      }
      ++idxCol;
    }

    // weighting with sample size
    // what rf5 of breiman/cutler does
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // calc. the true dec. of gini.
    //decOfGini = nodeGini - 1 + decOfGini / ntwmofbest;


    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

  template<class T>
  static ClassAtom<T, uli_t> *CARTginiSrtWithMissings(
      FittingFctPar<T> &fitFctPar) {
    INode<T> *&parent = fitFctPar.parent;
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    uli_t iteration = fitFctPar.iteration;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxMissVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    T cutoffT, cutoffTlast;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini;
    double ntwm, ntwmofbest;

    a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;
    k = j = 0;
    cutoffT = cutoffTlast = 0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);

    // init count
    idxRow = rowMaskVec.begin();
    while (idxRow != rowMaskVec.end()) {
      nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    Nt = bd;
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
#endif
      gini = 0.0;

      // get column data information
      data.getRowMaskedColNoMissing2(
          *idxCol, rowMaskVec, colVec, idxVec, idxMissVec);

      // if data contains only missings then skip the current column
      if (colVec.size() == 0) {
        ++idxCol;
        continue;
      }

      // decount missings
      nr = nrinit;
      ntwmofbest = ntwm = Nt;
      idxRow = idxMissVec.begin();
      while (idxRow != idxMissVec.end()) {
        nr[data.depIdxVec[*idxRow]] -= rowWeight[*idxRow];
        ntwm -= rowWeight[*idxRow];
        ++idxRow;
      }

      Helper::makePairVec<T, uli_t>(colVec, idxVec, rowPairVec);
      sort(rowPairVec.begin(), rowPairVec.end());

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }

      // iterate over all cutoffs in rowPairVec
      cutoffT = missingcode;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        // change cutoff if necessary
        // and calc gini index
        if (cutoffT != itRowPairVec->first) {
          // first iteration check
          (missingcode == cutoffT) ? cutoffTlast = itRowPairVec->first
              : cutoffTlast = cutoffT;
          cutoffT = itRowPairVec->first;

          // calculate updated gini
          // left child.
          itlr = nl.begin();
          a = 0.0;
          ad = 0.0;
          while (itlr != nl.end()) {
            a += (*itlr) * (*itlr);
            ad += (*itlr);
            ++itlr;
          }
          (ad != 0.0) ? (a /= ad)
              : (a = 0.0);

          // right child.
          itlr = nr.begin();
          b = 0.0;
          bd = 0.0;
          while (itlr != nr.end()) {
            b += (*itlr) * (*itlr);
            bd += (*itlr);
            ++itlr;
          }
          (bd != 0.0) ? (b /= bd)
              : (b = 0.0);

          gini = a + b;

#ifdef FITTINGFCTDEBUG2
          Helper::printVec<double>(nl);
          Helper::printVec<double>(nr);
          Helper::printVec<double>(rowWeight);
          Helper::printVec<T>(colVec);
          Helper::printVec<T>(classVarVec);
          Helper::printVec<uli_t>(rowMaskVec);
          *fitFctPar.io->outVerbose << (double) data.at(0, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(3, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(4, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << (double) data.at(8, *idxCol) << std::endl;
          *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
          *fitFctPar.io->outVerbose << "gini " << gini << std::endl;
          *fitFctPar.io->outVerbose << "cutoffT " << cutoffT << std::endl;
          *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
          *fitFctPar.io->outVerbose << "giniBestIdx " << giniBestIdx << std::endl;
          *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
          *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif

          // is it a good gini check
					if ((gini > decOfGini) ||
							((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
            decOfGini = gini;
            giniBestIdx = *idxCol;
            bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
            ntwmofbest = ntwm;
          }
        }

#ifdef FITTINGFCTDEBUG
        *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
        *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
        *fitFctPar.io->outVerbose << "nodeGini " << nodeGini << std::endl;
        *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
        *fitFctPar.io->outVerbose << "(T)bestCutoff " << (T) bestCutoff << std::endl;
#endif

        // get the differences
        k = data.depIdxVec[itRowPairVec->second];
        nl[k] += rowWeight[itRowPairVec->second];
        nr[k] -= rowWeight[itRowPairVec->second];

        ++itRowPairVec;
      }
      ++idxCol;
    }

    // weighting with sample size
    // what rf5 of breiman/cutler do
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // calc. the true dec. of gini.
    //decOfGini = nodeGini - 1 + decOfGini / ntwmofbest;

    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    fitFctPar.importance->add(classifier->getVarID(), decOfGini);

    return classifier;
  }

  /*
   * Trend-Test-CART
   *
   */
  template<class T>
  static ClassAtom<T, uli_t> *trendCARTgini(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &bestIdx = *fitFctPar.bestVar;

    uli_t j, k, giniBestIdx;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<T> colVec, classVarVec, giniVec, outcomeVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxMissVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff = 0;
    double a, b, bno, bde, ad, bd, Nt;
    double gini, stdDecOfGini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double ntwm, ntwmofbest;
    T cutoffT, cutoffTlast;

    if (rowMaskVec.size() <= 2) {
      decOfGini = 0;
      return classifier;
    }

    giniBestIdx = 0;
    a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;
    k = j = 0;
    cutoffT = cutoffTlast = 0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();

    // 2 x 5 conting. table (see LOTUS paper by Loh)
    // (aff. and unaff.) x (counts of X-values grouped in 5 sets)
    uli_t J = 5;
    std::vector<T> quantVec;

    // quantile indexer
    uli_t i;
    typename std::vector<T>::iterator itq;
    std::vector<double> nReset(J, 0);
    std::vector<double> nAff;
    std::vector<double> nUnaff;
    std::vector<T> xj;
    double n = 0;
    double XL, bestXL;
    double xweighted;

    bestXL = 0.0;

    // the name "aff." or "unaff." must not accord to the real aff. status
    T afftype = data.at(rowMaskVec[0], depVar);

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      gini = 0.0;

      //*fitFctPar.io->outVerbose << *idxCol << "," << std::endl;

      // get variable's quantiles
      data.getRowMaskedColNoMissing(*idxCol, rowMaskVec, colVec, idxVec);

      // if data contains only missings then skip the current column
      if (colVec.size() == 0) {
        ++idxCol;
        continue;
      }

      Helper::makePairVec(colVec, idxVec, rowPairVec);
      sort(rowPairVec.begin(), rowPairVec.end());

      // count samples
      n = 0;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        n += rowWeight[itRowPairVec->second];
        ++itRowPairVec;
      }

      // get quantiles
      quantVec.clear();
      quantVec.push_back(0);
      i = 0;
      j = 0;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        if (i > (uli_t) floor(n * (j + 1) / J)) {
          //quantVec[j] = colVec[i];
          quantVec[j] = itRowPairVec->first;
          quantVec.push_back(0);
          ++j;
        }
        i += (uli_t) rowWeight[itRowPairVec->second];
        ++itRowPairVec;
      }

      quantVec[j] = *(colVec.rbegin());

      // calc trend test score
      xj = quantVec;
      nAff = nUnaff = nReset;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        for (i = 0; i < quantVec.size(); ++i) {
          if (i == 0) {
            if (data.at(itRowPairVec->second, *idxCol) <= quantVec[0]) {
              if (data.at(itRowPairVec->second, depVar) == afftype)
                nAff[0] += rowWeight[itRowPairVec->second];
              else
                nUnaff[0] += rowWeight[itRowPairVec->second];
              break;
            }
          } else {
            if (data.at(itRowPairVec->second, *idxCol) > quantVec[i - 1]
                && data.at(itRowPairVec->second, *idxCol) <= quantVec[i]) {
              if (data.at(itRowPairVec->second, depVar) == afftype)
                nAff[i] += rowWeight[itRowPairVec->second];
              else
                nUnaff[i] += rowWeight[itRowPairVec->second];
              break;
            }
          }
        }
        ++itRowPairVec;
      }

      xweighted = 0.0;
      for (j = 0; j < J; ++j) {
        xweighted += xj[j] * (nAff[j] + nUnaff[j]);
      }

      bno = 0.0;
      for (j = 0; j < J; ++j) {
        bno += (nAff[j] + nUnaff[j]) * (nUnaff[j] - Helper::getSum<double>(
            nUnaff)) / n * (xj[j] - xweighted / n);
      }

      bde = 0.0;
      for (j = 0; j < J; ++j) {
        bde += (nAff[j] + nUnaff[j]) * (xj[j] - xweighted / n) * (xj[j]
            - xweighted / n);
      }

      /*
       if (bde == 0.0) {
       //++idxCol;
       //continue;
       b = 0.0;
       } else {
       b = bno / bde;
       }
       */

      XL = 0.0;
      for (j = 0; j < J; ++j) {
        bno = nUnaff[j] / n - (Helper::getSum<double>(nUnaff) / n + b * (xj[j]
            - xweighted / n));
        XL += (nAff[j] + nUnaff[j]) * bno * bno * n * n
            / Helper::getSum<double>(nAff) / Helper::getSum<double>(nUnaff);
      }

      if (XL == XL) {
				if ((XL > bestXL) ||
						((gsl_rng_uniform_int(data.par.rng, 2)==0) && (XL == bestXL))){
          bestXL = XL;
          bestIdx = *idxCol;
          bestCutoff = xj[(J - 1) / 2];
        }
      }
      ++idxCol;

    }

    // bonf. cor. p-value
    double pval = gsl_cdf_chisq_Q(bestXL, J - 2) * colMaskVec.size();

    // stop if diff. of best variable is not signif.
    if ((fitFctPar.depth != 0) && (pval > 0.05)) {
      decOfGini = 0.0;
      return NULL;
    }

    /*
     if (parentsBest == bestIdx) {
     decOfPurity = 0.0;
     return NULL;
     }
     */
    // get best split point in variable
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);

    // init count
    idxRow = rowMaskVec.begin();
    while (idxRow != rowMaskVec.end()) {
      nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    Nt = bd;
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;

    //get best column
    data.getRowMaskedColNoMissing2(
        bestIdx, rowMaskVec, colVec, idxVec, idxMissVec);

    // if data contains only missings then skip the current column
    if (colVec.size() == 0) {
      decOfGini = 0;
      classifier->setVarID(bestIdx);
      classifier->setThreshold(MISSINGCODE);
      classifier->setMissingCode(data.getMissingCode());
      return classifier;
    }

    // decount missings
    nr = nrinit;
    ntwmofbest = ntwm = Nt;
    idxRow = idxMissVec.begin();
    while (idxRow != idxMissVec.end()) {
      nr[data.depIdxVec[*idxRow]] -= rowWeight[*idxRow];
      ntwm -= rowWeight[*idxRow];
      ++idxRow;
    }

    Helper::makePairVec<T, uli_t>(colVec, idxVec, rowPairVec);
    sort(rowPairVec.begin(), rowPairVec.end());

    // initial reset of the left node
    itlr = nl.begin();
    while (itlr != nl.end()) {
      *itlr = 0.0;
      ++itlr;
    }

    // iterate over all cutoffs in rowPairVec
    cutoffT = missingcode;
    itRowPairVec = rowPairVec.begin();
    while (itRowPairVec != rowPairVec.end()) {
      // change cutoff if necessary
      // and calc gini index
      if (cutoffT != itRowPairVec->first) {
        // first iteration check
        (missingcode == cutoffT) ? cutoffTlast = itRowPairVec->first
            : cutoffTlast = cutoffT;
        cutoffT = itRowPairVec->first;

        // calculate updated gini
        // left child.
        itlr = nl.begin();
        a = 0.0;
        ad = 0.0;
        while (itlr != nl.end()) {
          a += (*itlr) * (*itlr);
          ad += (*itlr);
          ++itlr;
        }
        (ad != 0.0) ? (a /= ad)
            : (a = 0.0);

        // right child.
        itlr = nr.begin();
        b = 0.0;
        bd = 0.0;
        while (itlr != nr.end()) {
          b += (*itlr) * (*itlr);
          bd += (*itlr);
          ++itlr;
        }
        (bd != 0.0) ? (b /= bd)
            : (b = 0.0);

        gini = a + b;

        // is it a good gini check
        if ((gini > decOfGini) ||
						((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
          decOfGini = gini;
          bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
          ntwmofbest = ntwm;
        }
      }

      // get the differences
      k = data.depIdxVec[itRowPairVec->second];
      nl[k] += rowWeight[itRowPairVec->second];
      nr[k] -= rowWeight[itRowPairVec->second];

      ++itRowPairVec;
    }

    // output dec. of gini
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // build classifier
    classifier->setVarID(bestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

  /*
   * LOTUS
   *
   * CAUTION:
   * This fct expects only two classes in the dependent variable (Y).
   */
  template<class T>
  static ClassAtom<T, uli_t> *LOTUS(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    double nodeDev = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfPurity = *fitFctPar.treeImprovement;
    uli_t &bestIdx = *fitFctPar.bestVar;

    uli_t j;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    typename std::vector<std::pair<T, uli_t> >::reverse_iterator itRevRowPairVec;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<double> nl, nr, nrinit;
		std::vector<T > leftNodeVec, rightNodeVec, leftOutcomeVec, rightOutcomeVec;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier;
    classifier = new T2ClassAtom<T> ();
    double bestCutoff = 0;
    double a, b, bno, bde, ad, bd;
    uli_t parentsBest = bestIdx;
    typename std::vector<T>::iterator itT;
    typename std::vector<T>::reverse_iterator itRevT;

    a = b = ad = bd = 0.0;
    //		decOfGini = -std::numeric_limits<double>::infinity();

    // 2 x 5 conting. table (see LOTUS paper by Loh)
    // (aff. and unaff.) x (counts of X-values grouped in 5 sets)
    uli_t J = 5;
    std::vector<T> quantVec;

    // quantile indexer
    uli_t i;
    typename std::vector<T>::iterator itq;
    std::vector<double> nReset(J, 0);
    std::vector<double> nAff;
    std::vector<double> nUnaff;
    std::vector<T> xj;
    double n = 0;
    double XL, bestXL;
    double xweighted;

    bestXL = -1.0;

    // the name "aff." or "unaff." must not accord to the real aff. status
    T afftype = data.at(rowMaskVec[0], depVar);

		// 0) search best variable
    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      // get variable's quantiles
      data.getRowMaskedColNoMissing(*idxCol, rowMaskVec, colVec, idxVec);

      // if data contains only missings then skip the current column
      if (colVec.size() == 0) {
        ++idxCol;
        continue;
      }

      Helper::makePairVec(colVec, idxVec, rowPairVec);
      sort(rowPairVec.begin(), rowPairVec.end());

      // count samples
      n = 0;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
				//n += rowWeight[itRowPairVec->second];
        ++n;
        ++itRowPairVec;
      }

      // get quantiles
      quantVec.clear();
      quantVec.push_back(0);
      i = 0;
      j = 0;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        //if (i > (uli_t) floor(n * (j + 1) / J)) {
        if (i > (uli_t) floor(rowPairVec.size() * (j + 1) / J)) {
          //quantVec[j] = colVec[i];
          quantVec[j] = itRowPairVec->first;
          quantVec.push_back(0);
          ++j;
        }
        //i += (uli_t) rowWeight[itRowPairVec->second];
        ++i;
        ++itRowPairVec;
      }

      quantVec[j] = *(colVec.rbegin());

      // calc trend test score
      xj = quantVec;
      nAff = nUnaff = nReset;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        for (i = 0; i < quantVec.size(); ++i) {
          if (i == 0) {
            if (data.at(itRowPairVec->second, *idxCol) <= quantVec[0]) {
              if (data.at(itRowPairVec->second, depVar) == afftype)
                //nAff[0] += rowWeight[itRowPairVec->second];
                ++nAff[0];
              else
                //nUnaff[0] += rowWeight[itRowPairVec->second];
                ++nUnaff[0];
              break;
            }
          } else {
            if (data.at(itRowPairVec->second, *idxCol) > quantVec[i - 1]
                && data.at(itRowPairVec->second, *idxCol) <= quantVec[i]) {
              if (data.at(itRowPairVec->second, depVar) == afftype)
                //nAff[i] += rowWeight[itRowPairVec->second];
                ++nAff[i];
              else
                //nUnaff[i] += rowWeight[itRowPairVec->second];
                ++nUnaff[i];
              break;
            }
          }
        }
        ++itRowPairVec;
      }

      xweighted = 0.0;
      for (j = 0; j < J; ++j) {
        xweighted += xj[j] * (nAff[j] + nUnaff[j]);
      }

      bno = 0.0;
      for (j = 0; j < J; ++j) {
        bno += (nAff[j] + nUnaff[j]) * (nUnaff[j] - Helper::getSum<double>(
            nUnaff)) / n * (xj[j] - xweighted / n);
      }

      bde = 0.0;
      for (j = 0; j < J; ++j) {
        bde += (nAff[j] + nUnaff[j]) * (xj[j] - xweighted / n) * (xj[j]
            - xweighted / n);
      }

      b = bno / bde;

      XL = 0.0;
      for (j = 0; j < J; ++j) {
        bno = nUnaff[j] / n - (Helper::getSum<double>(nUnaff) / n + b * (xj[j]
            - xweighted / n));
        XL += (nAff[j] + nUnaff[j]) * bno * bno * n * n
            / Helper::getSum<double>(nAff) / Helper::getSum<double>(nUnaff);
      }

			if ((XL > bestXL) ||
					((gsl_rng_uniform_int(data.par.rng, 2)==0) && (XL == bestXL))){
        bestXL = XL;
        bestIdx = *idxCol;
        bestCutoff = xj[(J - 1) / 2];
      }
      ++idxCol;

    }

    if (parentsBest == bestIdx) {
      decOfPurity = 0.0;
      return NULL;
    }


		// 1) search split point
		
		// get quantiles
		data.getRowMaskedColNoMissing2AndPair(bestIdx, rowMaskVec, rowPairVec);
		sort(rowPairVec.begin(), rowPairVec.end());
		Helper::getQuantilesLotusSplit(rowPairVec, quantVec);
		quantVec.resize(unique(quantVec.begin(), quantVec.end()) - quantVec.begin());

		// prepare for split search
		Helper::putPairFirstToVec(rowPairVec, leftNodeVec);
		data.getRowMaskedColNoMissing(depVar, rowMaskVec, leftOutcomeVec);

		// 1.1) foreach 0.3, 0.4, 0.5, 0.6, and 0.7 quantiles of X
		//FILE *outputFile = fopen("stuff", "w");
		//fclose(outputFile);
		double dev, devBest;
		T cutoffBest;
		dev = i = cutoffBest = 0;
		devBest = std::numeric_limits<double>::infinity();
		itRevRowPairVec = rowPairVec.rbegin();
		itRevT = quantVec.rbegin();

		while (itRevRowPairVec != rowPairVec.rend()) {
			
			// next split/quantile? then ->
			// LOTUS(S)                                
			// a best simple linear regression model to every node
			// no multiple regression, because of huge data
			if (itRevRowPairVec->first <= *itRevT) {

				//   1.1.2) perform lr on left data and save deviance
				//          if it doesnt converged then jump to step 2)
				dev = Helper::performLrAndGetDeviance(leftNodeVec, leftOutcomeVec);

				//   1.1.3) perform lr on right data and save deviance
				//          if it doesnt converged then jump to step 2)
				dev += Helper::performLrAndGetDeviance(rightNodeVec, rightOutcomeVec);

				//   1.1.4) add both left and right deviance
				//          if it is the smallest til now
				//          then store it as the best
				if ((dev > devBest) ||
						((gsl_rng_uniform_int(data.par.rng, 2)==0) && (dev == devBest))){
					devBest = dev;
					cutoffBest = *itRevT;
				}
        
				++itRevT;
				if (itRevT == quantVec.rend()) break;
			}


			// put one sample from left to right node
			rightNodeVec.push_back(leftNodeVec.back());
			leftNodeVec.pop_back();
			rightOutcomeVec.push_back(leftOutcomeVec.back());
			leftOutcomeVec.pop_back();

			++itRevRowPairVec;
		}

		// 2) if dev==inf then perform lr to current node and convert it to terminal node
		if (devBest == std::numeric_limits<double>::infinity()) {
      fitFctPar.stopGrowing = true;
      return classifier;
		}

		// 3) return split point etc and leave the fitting function

    // build classifier
    classifier->setVarID(bestIdx);
    classifier->setThreshold(cutoffBest);
    classifier->setMissingCode(data.getMissingCode());
		
    //add importance
    decOfPurity = nodeDev - devBest; 
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfPurity / data.par.ntree);

    return classifier;
  }

		/*
		 * check data, if LOTUSginiHybrid is usable
		 */
		template<class T>
			static void checkUsabilityLOTUS(DataFrame<T> &data) {
				uli_t depVar = data.par.depVar;

    if (data.par.memMode > 1) {
      throw Exception(ERRORCODE_64);
    }

    if (data.varClassIndexer.size() <= depVar) {
      throw Exception(ERRORCODE_19);
    }

    if ((data.varClassIndexer[depVar][0] == 0
        && data.varClassIndexer[depVar][1] == 1)
        || (data.varClassIndexer[depVar][0] == 1
            && data.varClassIndexer[depVar][1] == 0)) {
    } else {
      throw Exception(ERRORCODE_20);
    }
  }

  template<class T>
  static ClassAtom<T, uli_t> *SCARTgini(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    typename std::vector<T>::iterator itT;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    std::vector<uli_t>::iterator idxCol;
    TMClassAtom<T> *classifier = new TMClassAtom<T> ();
    //double speed = 0.01;

    decOfGini = std::numeric_limits<double>::infinity();
    double gini;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      data.getRowMaskedColNoMissing(*idxCol, rowMaskVec, colVec, maskVec);
      //categorical
      // search lowest gini for each threshold

      gini = 0.0;
      for (j = 0; j < classVarVec.size(); ++j) {
        Helper::getMask<T>(classVarVec[j], colVec, maskVec);
        Helper::getIndexedVector<uli_t>(rowMaskVec, maskVec, rowMaskForDepVec);
        gini += (double) maskVec.size() * Helper::purityGini<T>(
            data, depVar, rowMaskForDepVec, rowWeight);

      }

      gini /= colVec.size();

			if ((gini > decOfGini) ||
					((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
				decOfGini = gini;
				giniBestIdx = *idxCol;
      }

      ++idxCol;
    }

    decOfGini = nodeGini - decOfGini;

    // build classifier
    data.getRowMaskedColNoMissing(giniBestIdx, rowMaskVec, colVec, maskVec);
    Helper::getClasses<T>(colVec, classVarVec);
    sort(classVarVec.begin(), classVarVec.end());

    if (classVarVec.size() > 1) {
      itT = classVarVec.begin();
      ++itT;
      while (itT != classVarVec.end()) {
        *(itT - 1) = (T) ((*itT + *(itT - 1)) / 2.0);
        ++itT;
      }
      classVarVec.pop_back();
    }

    classifier->setVarID(giniBestIdx);
    classifier->setThresholdsNoSort(classVarVec);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

  /*
   * unbiasedCARTgini
   *
   * Find the best numerical/ordinal variable with categorical outcome.
   * The gini method is described in CART (Leo Breiman, et al. 1983).
   * It regards missing values.
   * unbiased extension by Daniel Schwarz
   */
  template<class T>
  static ClassAtom<T, uli_t> *unbiasedCARTgini__(FittingFctPar<T> &fitFctPar) {
    INode<T> *&parent = fitFctPar.parent;
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    uli_t iteration = fitFctPar.iteration;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> dumpVec;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff = 0;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    double cutoff;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini = 0.0;
    double ntwm, ntwmofbest;

    a = b = ad = bd = 0.0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();

    // 2 x 5 conting. table
    // (aff. and unaff.) x (counts of X-values grouped in 5 sets)
    uli_t J = 2;
    std::vector<T> quantVec;
    uli_t i;
    typename std::vector<T>::iterator itq;
    std::vector<double> nReset(J, 0);
    std::vector<double> nAff;
    std::vector<double> nUnaff;
    double n = 0;

    // initial count (regard weights from bootstrap)
    k = j = 0;
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);
    idxVec.clear();
    idxRow = rowMaskVec.begin();
    j = 0;
    Nt = 0;
    while (idxRow != rowMaskVec.end()) {
      idxVec.push_back(data.depIdxVec[*idxRow]);
      nrinit[idxVec[j]] += rowWeight[*idxRow];
      Nt += rowWeight[*idxRow];
      ++j;
      ++idxRow;
    }
    ntwmofbest = Nt;

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    b = 0.0;
    bd = 0.0;
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;
    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

      gini = 0.0;

      // get variable's classes
      /*
       data.getRowMaskedCol(*idxCol, rowMaskVec, colVec);
       Helper::getClasses<T >(colVec, classVarVec);
       sort(classVarVec.begin(), classVarVec.end());

       */

      // get variable's quantiles
      data.getRowMaskedColNoMissing(*idxCol, rowMaskVec, colVec, dumpVec);

      // if data contains only missings then skip the current column
      if (colVec.size() == 0) {
        ++idxCol;
        continue;
      }

      Helper::makePairVec<T, uli_t>(colVec, dumpVec, rowPairVec);
      sort(rowPairVec.begin(), rowPairVec.end());

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "colVec ";
      for (i = 0; i < colVec.size(); i++) {
        *fitFctPar.io->outVerbose << colVec[i] << " ";
      }
      *fitFctPar.io->outVerbose << std::endl;
#endif
      // count samples
      n = 0;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        n += rowWeight[itRowPairVec->second];
        ++itRowPairVec;
      }

      // get quantiles
      quantVec.clear();
      quantVec.push_back(0);
      i = 0;
      quantVec[0] = rowPairVec.begin()->first;
      quantVec.push_back(0);
      j = 1;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        if (i > (uli_t) floor(n * j / (J - 1))) {
          //quantVec[j] = colVec[i];
          quantVec[j] = itRowPairVec->first;
          quantVec.push_back(0);
          ++j;
        }
        i += (uli_t) rowWeight[itRowPairVec->second];
        ++itRowPairVec;
      }
      quantVec[j] = rowPairVec.rbegin()->first;

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "j " << j << std::endl;
      *fitFctPar.io->outVerbose << "quantVec.size() " << quantVec.size() << std::endl;
      *fitFctPar.io->outVerbose << "quantVec ";
      for (i = 0; i < quantVec.size(); i++) {
        *fitFctPar.io->outVerbose << quantVec[i] << " ";
      }
      *fitFctPar.io->outVerbose << std::endl;
#endif

      Helper::getClasses<T>(quantVec, classVarVec);
      sort(classVarVec.begin(), classVarVec.end());

      // break if there is just one class
      if (classVarVec.size() == 1 && classVarVec[0] == missingcode) {
        ++idxCol;
        continue;
      } else if ((classVarVec.size() <= 2 && classVarVec[0] == missingcode)
          || (classVarVec.size() <= 2 && classVarVec[1] == missingcode)
          || classVarVec.size() == 1) {
        gini = stdDecOfGini;
        if (gini > decOfGini) {
          decOfGini = gini;
          giniBestIdx = *idxCol;
          if (classVarVec.size() == 1 || (classVarVec.size() <= 2
              && classVarVec[1]) == missingcode)
            bestCutoff = (double) classVarVec[0];
          else
            bestCutoff = (double) classVarVec[1];
        }
        ++idxCol;
        continue;
      }

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }
      nr = nrinit;
      ntwm = Nt; //"degrade" missing samples

      // iterate over all cutoffs
      for (j = 0; j < classVarVec.size() - 1; ++j) {
        // intrinsic cutoff between two value which is
        // (v_i + v_{i+1})/2
        cutoff = (double) (classVarVec[j]);

        k = 0;
        idxRow = rowMaskVec.begin();
        if (classVarVec[j] == missingcode) {
          // degrade missings
          while (idxRow != rowMaskVec.end()) {
            if (data.at(*idxRow, *idxCol) == missingcode) {
              ntwm -= rowWeight[*idxRow];
            }
            ++k;
            ++idxRow;
          }
          continue;
        } else {
          // get the differences
          while (idxRow != rowMaskVec.end()) {
            if (data.at(*idxRow, *idxCol) == cutoff) {
              nl[idxVec[k]] += rowWeight[*idxRow];
              nr[idxVec[k]] -= rowWeight[*idxRow];
            }
            ++k;
            ++idxRow;
          }
        }

        // calculate updated gini
        // left child.
        itlr = nl.begin();
        a = 0.0;
        ad = 0.0;
        while (itlr != nl.end()) {
          a += (*itlr) * (*itlr);
          ad += (*itlr);
          ++itlr;
        }
        (ad != 0.0) ? (a /= ad)
            : (a = 0.0);

        // right child.
        itlr = nr.begin();
        b = 0.0;
        bd = 0.0;
        while (itlr != nr.end()) {
          b += (*itlr) * (*itlr);
          bd += (*itlr);
          ++itlr;
        }
        (bd != 0.0) ? (b /= bd)
            : (b = 0.0);

        gini = a + b;

        // is it a good gini?
				if ((gini > decOfGini) ||
						((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
					decOfGini = gini;
          giniBestIdx = *idxCol;
          bestCutoff = (cutoff + (double) classVarVec[j + 1]) / 2;
          ntwmofbest = ntwm;
        }
      }
      ++idxCol;
    }

    // output dec. of gini
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    fitFctPar.importance->add(classifier->getVarID(), decOfGini);

#ifdef FITTINGFCTDEBUG

    *fitFctPar.io->outVerbose << "> decOfGini " << decOfGini << std::endl;
    *fitFctPar.io->outVerbose << "nodeGini: " << nodeGini << std::endl;
    *fitFctPar.io->outVerbose << "bestCutoff: " << bestCutoff << std::endl;
    *fitFctPar.io->outVerbose << "ntwmofbest: " << ntwmofbest << std::endl;
    *fitFctPar.io->outVerbose << "quantVec: " << &quantVec << std::endl;
    *fitFctPar.io->outVerbose << "classifier: " << classifier << std::endl;
    *fitFctPar.io->outVerbose << "sizeof(classifier): " << sizeof(ClassAtom<T, uli_t> )
    << std::endl;
    if (classifier != NULL)
    *fitFctPar.io->outVerbose << "classifier: " << *classifier << std::endl;
    *fitFctPar.io->outVerbose << "ntwmofbest * (nodeGini - 1): " << ntwmofbest * (nodeGini - 1)
    << std::endl;
#endif

    return classifier;
  }

  template<class T>
  static ClassAtom<T, uli_t> *unbiasedCARTgini(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t i, j, k, J;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<T> quantVec;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxMissVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    T cutoffT, cutoffTlast;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini;
    double ntwm, ntwmofbest;

    a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;
    i = j = k = 0;
    cutoffT = cutoffTlast = 0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);
    // parts
    J = 3;

    // init count
    idxRow = rowMaskVec.begin();
    while (idxRow != rowMaskVec.end()) {
      nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    Nt = bd;
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
#endif
      gini = 0.0;

      // get column data information
      data.getRowMaskedColNoMissing2(
          *idxCol, rowMaskVec, colVec, idxVec, idxMissVec);

      // if data contains only missings then skip the current column
      if (colVec.size() == 0) {
        ++idxCol;
        continue;
      }

      // decount missings
      nr = nrinit;
      ntwmofbest = ntwm = Nt;
      idxRow = idxMissVec.begin();
      while (idxRow != idxMissVec.end()) {
        nr[data.depIdxVec[*idxRow]] -= rowWeight[*idxRow];
        ntwm -= rowWeight[*idxRow];
        ++idxRow;
      }

      // sort variable content
      Helper::makePairVec<T, uli_t>(colVec, idxVec, rowPairVec);
      sort(rowPairVec.begin(), rowPairVec.end());

      /* get quantiles */
      // init
      quantVec.clear();
      quantVec.push_back(0);
      quantVec[0] = rowPairVec.begin()->first;
      quantVec.push_back(0);

      // set quantVec
      i = 0;
      j = 1;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        if (i > (uli_t) floor(ntwm * j / (J - 1))) {
          quantVec[j] = itRowPairVec->first;
          quantVec.push_back(0);
          ++j;
        }
        i += (uli_t) rowWeight[itRowPairVec->second];
        ++itRowPairVec;
      }
      quantVec[j] = rowPairVec.rbegin()->first;

      // filter and sort quantiles
      Helper::getClasses<T>(quantVec, classVarVec);
      sort(classVarVec.begin(), classVarVec.end());

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }

      // iterate over all cutoffs in rowPairVec
      j = 0;
      cutoffT = missingcode;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        // change cutoff if necessary and calc gini index
        if (missingcode == cutoffT || cutoffT < itRowPairVec->first) {
          // first iteration check
          (missingcode == cutoffT) ? cutoffTlast = classVarVec[j]
              : cutoffTlast = cutoffT;
          cutoffT = classVarVec[j];
          ++j;

          // calculate updated gini
          // left child.
          itlr = nl.begin();
          a = 0.0;
          ad = 0.0;
          while (itlr != nl.end()) {
            a += (*itlr) * (*itlr);
            ad += (*itlr);
            ++itlr;
          }
          (ad != 0.0) ? (a /= ad)
              : (a = 0.0);

          // right child.
          itlr = nr.begin();
          b = 0.0;
          bd = 0.0;
          while (itlr != nr.end()) {
            b += (*itlr) * (*itlr);
            bd += (*itlr);
            ++itlr;
          }
          (bd != 0.0) ? (b /= bd)
              : (b = 0.0);

          gini = a + b;

          // is it a good gini check
					if ((gini > decOfGini) ||
							((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
						decOfGini = gini;
            giniBestIdx = *idxCol;
            bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
            ntwmofbest = ntwm;
          }
        }

        // get the differences
        k = data.depIdxVec[itRowPairVec->second];
        nl[k] += rowWeight[itRowPairVec->second];
        nr[k] -= rowWeight[itRowPairVec->second];

        ++itRowPairVec;
      }
      ++idxCol;
    }

    // weighting with sample size
    // what rf5 of breiman/cutler do
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // calc. the true dec. of gini.
    //decOfGini = nodeGini - 1 + decOfGini / ntwmofbest;

    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

  /*
   * CARTginiSrt
   *
   * Find the best numerical/ordinal variable with categorical outcome.
   * The gini method is described in CART (Leo Breiman, et al. 1983).
   * It regards missing values.
   * Idea of original algorithm from A. Cutler
   */
  template<class T>
  static ClassAtom<T, uli_t> *imputeTree(FittingFctPar<T> &fitFctPar) {
    INode<T> *&parent = fitFctPar.parent;
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    uli_t iteration = fitFctPar.iteration;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<std::pair<T, uli_t> > rowPairVec;
    typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> idxMissVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    T cutoffT, cutoffTlast;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini;
    double ntwm, ntwmofbest;

    a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;
    k = j = 0;
    cutoffT = cutoffTlast = 0;
    decOfGini = -std::numeric_limits<double>::infinity();
    T missingcode = data.getMissingCode();
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);

    // init count
    idxRow = rowMaskVec.begin();
    while (idxRow != rowMaskVec.end()) {
      nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    Nt = bd;
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;

    idxCol = colMaskVec.begin();
    while (idxCol != colMaskVec.end()) {
      //skip dep var
      if (*idxCol == depVar) {
        ++idxCol;
        continue;
      }

#ifdef FITTINGFCTDEBUG
      *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
#endif
      gini = 0.0;

      // get column data information
      data.getRowMaskedColNoMissing2(
          *idxCol, rowMaskVec, colVec, idxVec, idxMissVec);

      // if data contains only missings then skip the current column
      if (colVec.size() == 0) {
        ++idxCol;
        continue;
      }

      // decount missings
      nr = nrinit;
      ntwmofbest = ntwm = Nt;
      idxRow = idxMissVec.begin();
      while (idxRow != idxMissVec.end()) {
        nr[data.depIdxVec[*idxRow]] -= rowWeight[*idxRow];
        ntwm -= rowWeight[*idxRow];
        ++idxRow;
      }

      Helper::makePairVec<T, uli_t>(colVec, idxVec, rowPairVec);
      sort(rowPairVec.begin(), rowPairVec.end());

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }

      // iterate over all cutoffs in rowPairVec
      cutoffT = missingcode;
      itRowPairVec = rowPairVec.begin();
      while (itRowPairVec != rowPairVec.end()) {
        // change cutoff if necessary
        // and calc gini index
        if (cutoffT != itRowPairVec->first) {
          // first iteration check
          (missingcode == cutoffT) ? cutoffTlast = itRowPairVec->first
              : cutoffTlast = cutoffT;
          cutoffT = itRowPairVec->first;

          // calculate updated gini
          // left child.
          itlr = nl.begin();
          a = 0.0;
          ad = 0.0;
          while (itlr != nl.end()) {
            a += (*itlr) * (*itlr);
            ad += (*itlr);
            ++itlr;
          }
          (ad != 0.0) ? (a /= ad)
              : (a = 0.0);

          // right child.
          itlr = nr.begin();
          b = 0.0;
          bd = 0.0;
          while (itlr != nr.end()) {
            b += (*itlr) * (*itlr);
            bd += (*itlr);
            ++itlr;
          }
          (bd != 0.0) ? (b /= bd)
              : (b = 0.0);

          gini = a + b;

          // weight gini with distance to variable of parent node
          if (parent != NULL) {
            //TODO: geht nicht mehr, wg. varID
            //gini *= pow(1 - (fabs(*idxCol - parent->varID) / data.par.ncol), 6);

          }

          // is it a good gini check
					if ((gini > decOfGini) ||
							((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
            decOfGini = gini;
            giniBestIdx = *idxCol;
            bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
            ntwmofbest = ntwm;
          }
        }

        // get the differences
        k = data.depIdxVec[itRowPairVec->second];
        nl[k] += rowWeight[itRowPairVec->second];
        nr[k] -= rowWeight[itRowPairVec->second];

        ++itRowPairVec;
      }
      ++idxCol;
    }
    // output dec. of gini
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    if (iteration == 10)
      decOfGini = 0;

    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T> *) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

  /*
   * CART gini. But first split is performed by variable in file "corrected.tree"
   */
  template<class T>
  static ClassAtom<T, uli_t> *CARTginiCor(FittingFctPar<T> &fitFctPar) {
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;
    std::vector<uli_t> &rowMaskVec = *fitFctPar.rowMaskVec;
    std::vector<uli_t> &colMaskVec = *fitFctPar.colMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    double nodeGini = fitFctPar.nodePerformance;
    //OUTPUT
    double &decOfGini = *fitFctPar.treeImprovement;
    uli_t &giniBestIdx = *fitFctPar.bestVar;

    uli_t j, k;
    std::vector<T> colVec;
    std::vector<T> classVarVec;
    std::vector<T> giniVec;
    std::vector<T> outcomeVec;
    std::vector<double> nl, nr, nrinit;
    std::vector<double>::iterator itlr;
    std::vector<uli_t> idxVec;
    std::vector<uli_t> maskVec;
    std::vector<uli_t> classMaskVec;
    std::vector<uli_t> rowMaskForDepVec;
    std::vector<uli_t> classOutcomeVec;
    typename std::vector<uli_t>::iterator idxRow;
    typename std::vector<uli_t>::iterator idxCol;
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestCutoff = 0;
    double a, b, ad, bd, Nt; // Nt = number of samples in current node
    T cutoff;
    double gini;
    uli_t depVarClassSize = data.varCategories[depVar].size();
    double stdDecOfGini = 0.0;
    double ntwm, ntwmofbest;

    /*
     * Assign plugin parameters
     */
    std::string pluginPar(data.par.pluginPar);

    /*
     * Optional feature: First node splits always variable xyz
     */
    bool isCorrectedTree = false;
    bool isFirstVar = (fitFctPar.depth == 0) ? true
        : false;
    uli_t idxFirstVar = 0;

    // Load 1st var.
    if ((isFirstVar) && (pluginPar != "")) {
      idxFirstVar = data.getIndexOfVarName(pluginPar);
      isCorrectedTree = true;
    }

    a = b = ad = bd = 0.0;
    decOfGini = -std::numeric_limits<double>::infinity();

    // initial count
    k = j = 0;
    nl.clear();
    nl.resize(depVarClassSize, 0.0);
    nrinit.clear();
    nrinit.resize(depVarClassSize, 0.0);
    idxVec.clear();
    idxRow = rowMaskVec.begin();
    j = 0;
    Nt = 0;
    while (idxRow != rowMaskVec.end()) {
      idxVec.push_back(data.depIdxVec[*idxRow]);
      nrinit[idxVec[j]] += rowWeight[*idxRow];
      ++j;
      ++idxRow;
    }

    // 1 + gini (without normalisation) of the whole node
    itlr = nrinit.begin();
    b = 0.0;
    bd = 0.0;
    while (itlr != nrinit.end()) {
      b += (*itlr) * (*itlr);
      bd += (*itlr);
      ++itlr;
    }
    (bd != 0.0) ? (b /= bd)
        : (b = 0.0);
    stdDecOfGini = b;
    Nt = bd;
    ntwmofbest = Nt;

    uli_t curVar;

    idxCol = colMaskVec.begin();
    while ((idxCol != colMaskVec.end()) || (isFirstVar)) {
      //skip dep var
      if (isFirstVar && isCorrectedTree) {
        curVar = idxFirstVar;
      } else {
        if ((*idxCol == depVar) || (*idxCol == idxFirstVar)) {
          ++idxCol;
          continue;
        }
        curVar = *idxCol;
      }

      gini = 0.0;

      // get variable's classes
      data.getRowMaskedCol(curVar, rowMaskVec, colVec);
      Helper::getClasses<T>(colVec, classVarVec);
      sort(classVarVec.begin(), classVarVec.end());
      //data.getRowMaskedClassesInCol(*idxCol, rowMaskVec, classVarVec);

      // break if there is just one class
      if (classVarVec.size() == 1) {
        gini = stdDecOfGini;
        if (gini > decOfGini) {
          decOfGini = gini;
          giniBestIdx = curVar;
          bestCutoff = (double) classVarVec[0];
        }
        if (isFirstVar) {
          break;
        } else {
          ++idxCol;
          continue;
        }
      }

      // initial reset of the left node
      itlr = nl.begin();
      while (itlr != nl.end()) {
        *itlr = 0.0;
        ++itlr;
      }
      nr = nrinit;
      ntwm = Nt; //"degrade" missing samples

      // iterate over all cutoffs
      for (j = 0; j < classVarVec.size() - 1; ++j) {
        // intrinsic cutoff between two value which is
        // (v_i + v_{i+1})/2
        cutoff = (classVarVec[j]);

        k = 0;
        idxRow = rowMaskVec.begin();
        // get the differences
        while (idxRow != rowMaskVec.end()) {
          if (data.at(*idxRow, curVar) == cutoff) {
            nl[idxVec[k]] += rowWeight[*idxRow];
            nr[idxVec[k]] -= rowWeight[*idxRow];
          }
          ++k;
          ++idxRow;
        }

        // calculate updated gini
        // left child.
        itlr = nl.begin();
        a = 0.0;
        ad = 0.0;
        while (itlr != nl.end()) {
          a += (*itlr) * (*itlr);
          ad += (*itlr);
          ++itlr;
        }
        (ad != 0.0) ? (a /= ad)
            : (a = 0.0);

        // right child.
        itlr = nr.begin();
        b = 0.0;
        bd = 0.0;
        while (itlr != nr.end()) {
          b += (*itlr) * (*itlr);
          bd += (*itlr);
          ++itlr;
        }
        (bd != 0.0) ? (b /= bd)
            : (b = 0.0);

        gini = a + b;

#ifdef FITTINGFCTDEBUG2
        Helper::printVec<double>(nl);
        Helper::printVec<double>(nr);
        Helper::printVec<double>(rowWeight);
        Helper::printVec<T>(colVec);
        Helper::printVec<T>(classVarVec);
        Helper::printVec<uli_t>(rowMaskVec);
        *fitFctPar.io->outVerbose << (double) data.at(0, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << (double) data.at(3, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << (double) data.at(4, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << (double) data.at(8, *idxCol) << std::endl;
        *fitFctPar.io->outVerbose << "*idxCol " << *idxCol << std::endl;
        *fitFctPar.io->outVerbose << "gini " << gini << std::endl;
        *fitFctPar.io->outVerbose << "cutoff " << cutoff << std::endl;
        *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
        *fitFctPar.io->outVerbose << "giniBestIdx " << giniBestIdx << std::endl;
        *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
        *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif
        // is it a good gini?
				if ((gini > decOfGini) ||
						((gsl_rng_uniform_int(data.par.rng, 2)==0) && (gini == decOfGini))){
          decOfGini = gini;
          giniBestIdx = curVar;
          bestCutoff = ((double) cutoff + (double) classVarVec[j + 1]) / 2.0;
          ntwmofbest = ntwm;
        }

      }
      if (isFirstVar) {
        break;
      } else {
        ++idxCol;
      }
    }
#ifdef FITTINGFCTDEBUG
    *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
    *fitFctPar.io->outVerbose << "ntwmofbest " << ntwmofbest << std::endl;
    *fitFctPar.io->outVerbose << "nodeGini " << nodeGini << std::endl;
    *fitFctPar.io->outVerbose << "bestCutoff " << bestCutoff << std::endl;
    *fitFctPar.io->outVerbose << "(T)bestCutoff " << (T) bestCutoff << std::endl;
#endif

    // weighting with sample size
    // what rf5 of breiman/cutler do
    decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

    // calc. the true dec. of gini.
    //decOfGini = nodeGini - 1 + decOfGini / ntwmofbest;

#ifdef FITTINGFCTDEBUG
    *fitFctPar.io->outVerbose << "decOfGini " << decOfGini << std::endl;
#endif
    // build classifier
    classifier->setVarID(giniBestIdx);
    classifier->setThreshold((T) bestCutoff);
    classifier->setMissingCode(data.getMissingCode());

    //add importance
    ((TImportance<T>*) fitFctPar.importance)->add(
        classifier->getVarID(), decOfGini / data.par.ntree);

    return classifier;
  }

};

#endif /*FITTINGFCT_H_*/
