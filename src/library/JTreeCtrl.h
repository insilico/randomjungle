#ifndef JTREECTRL_H_
#define JTREECTRL_H_

#ifdef __DEBUG__
#define JTREEDEBUG
#endif

/*
 * Includes
 */

#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "treedefs.h"
#include "Tree.h"
#include "INode.h"
#include "Exception.h"
#include "TermClassAtom.h"
#include "Importance.h"
#include "FittingFct.h"
#include "CondInf.h"
#include "DataTreeSet.h"
#include "RJunglePar.h"
#include "RJungleIO.h"
#include "RJungleGen.h"

#include "TimeProf.h"

/*
 * Source (Def. & Decl.)
 */

template<class T>
class JTreeCtrl {
public:
  JTreeCtrl() :
    data(NULL),
    colBuffer(NULL),
    rowBuffer(NULL),
    colChoose(NULL),
    rowChoose(NULL),
    importance(NULL)
  //      ,yaimp(NULL)
  {
  }
  ;

  JTreeCtrl(RJunglePar &par, RJungleIO &io, RJungleGen<T> &gen) :
    data(NULL),
    colBuffer(NULL),
    rowBuffer(NULL),
    colChoose(NULL),
    rowChoose(NULL),
    importance(NULL),
    //			yaimp(NULL),
        par(par), io(io), gen(gen) {
  }
  ;

  virtual ~JTreeCtrl() {
    if (colBuffer != NULL)
      delete[] colBuffer;
    if (rowBuffer != NULL)
      delete[] rowBuffer;
    if (importance != NULL)
      delete importance;
  }
  ;

  void addProximity2order(Node<T, uli_t> *node, DataFrame<double> *varProxiSum,
      std::vector<uli_t> &idxVec) {
    //    uli_t i, r, c, val;
    //
    //		for(i = 0; i < node->kids.size(); ++i) {
    //			if (i < node->kids.size() && node->kids.at(i) != NULL) {
    //        r = idxVec[node->varID];
    //        c = idxVec[node->kids.at(i)->varID];
    //        val = varProxiSum->at(r, c) + 1;
    //        varProxiSum->setArray(r, c, val);
    //  			addProximity2rder(
    //          node->kids.at(i),
    //          varProxiSum,
    //          idxVec
    //          );
    //			}
    //		}

  }

  virtual INode<T> *makeNode(DataFrame<T> &data, RJunglePar &par, std::vector<
      uli_t> &rowMaskVec, std::vector<double> &trainSetWeight, uli_t depth = 0,
      INode<T> *parent = NULL, double lastImprovement = 0) {

		TIMEPROF_START("JTreeCtrl~~makeNodeIntern");

    uli_t bestVar, i;
    double treeImprovement = -std::numeric_limits<double>::infinity();
    std::vector<T> colVec;
    std::vector<std::vector<uli_t> *> classMaskMultiVec;
    std::vector<uli_t> classVec, colMaskVec;
    ClassAtom<T, uli_t> *classifier;
    INode<T> *iNode = NULL;
    double improvementStop = 0.0; //10^-7
    double nodePerformance;
    i = bestVar = 0;

#ifdef JTREEDEBUG
    *io.outVerbose << "tree depth: " << depth << std::endl;
#endif

#ifdef __DEBUG2__
    *io.outVerbose << "data: " << std::endl << data << std::endl;
#endif

    //(parent != NULL) ? bestVar = parent->classigier->getVarID(): bestVar = 0;

    data.getRowMaskedCol(par.depVar, rowMaskVec, colVec);

    nodePerformance = gen.fct.getNodePerformance(data, par.depVar, rowMaskVec, trainSetWeight);

    /*
     if (depth > 15000) {
     //throw Exception(ERRORCODE_28);
     std::cerr << "Warning::28" << std::endl << std::flush;
     }
     */

#ifdef JTREEDEBUG
    *io.outVerbose << "node size: " << Helper::getSize<T>(colVec, data.getMissingCode()) << std::endl;
    *io.outVerbose << "sizeCutoff: " << this->sizeCutoff << std::endl;
    *io.outVerbose << "rowMaskVec: ";
    Helper::printVec(rowMaskVec);
    *io.outVerbose << "trainSetWeight: ";
    Helper::printVec(trainSetWeight);
    *io.outVerbose << "nodePerformance: " << nodePerformance << std::endl;
    *io.outVerbose << "par.maxTreeDepth: " << par.maxTreeDepth << std::endl;
#endif

    // randomly choose varPerNode variables per Node
    gsl_ran_choose(par.rng, colChoose, par.mtry, colBuffer, colBufferSize,
        sizeof(uli_t));

    gsl_ran_shuffle(par.rng, colChoose, par.mtry, sizeof(uli_t));

    for (uli_t m = 0; m < par.mtry; m++) {
      colMaskVec.push_back(colChoose[m]);
    }

    // break if nothing else to do
    // add terminal node
    if (Helper::getSize<T>(colVec, data.getMissingCode()) < this->sizeCutoff
        || nodePerformance == 0.0 || depth > par.maxTreeDepth) {

			TIMEPROF_STOP("JTreeCtrl~~makeNodeIntern");
      return new INode<T> (gen.fct.getTermResult(data, par, par.depVar,
						rowMaskVec, colMaskVec, trainSetWeight));
    }

    fitFctPar->parent = parent;
    fitFctPar->depth = depth;
    fitFctPar->rowMaskVec = &rowMaskVec;
    fitFctPar->colMaskVec = &colMaskVec;
    fitFctPar->nodePerformance = nodePerformance;
		fitFctPar->treeImprovement = &treeImprovement;
    fitFctPar->bestVar = &bestVar;
    fitFctPar->stopGrowing = false;
    fitFctPar->io = &(this->io);

		TIMEPROF_STOP("JTreeCtrl~~makeNodeIntern");
		TIMEPROF_START("JTreeCtrl~~FittingFct::getBestVariable");
    classifier = gen.fct.getBestVariable(*fitFctPar);
		TIMEPROF_STOP("JTreeCtrl~~FittingFct::getBestVariable");
		TIMEPROF_START("JTreeCtrl~~makeNodeIntern");

    //treeImprovement = fitFctPar.treeImprovement;
    //bestVar = fitFctPar.bestVar;

    // if decrease is not good enough stop building OR
    // if selected variable has got only one type of value
    if (
    //0.5 > gsl_rng_uniform(rng) ||
    fitFctPar->stopGrowing || fabs(treeImprovement) <= improvementStop
        || treeImprovement == std::numeric_limits<double>::infinity()
        || gen.fct.getNodePerformance(data, bestVar, rowMaskVec, trainSetWeight)
            == 0.0) {
			TIMEPROF_STOP("JTreeCtrl~~makeNodeIntern");
      return new INode<T>(gen.fct.getTermResult(data, par, par.depVar,
						rowMaskVec, colMaskVec, trainSetWeight));
    }

#ifdef JTREEDEBUG
    // save improvement results
    *io.outVerbose << "treeImprovement,";
    *io.outVerbose << "bestVar,";
    *io.outVerbose << "depth";
    *io.outVerbose << std::endl;
    *io.outVerbose << treeImprovement << ",";
    *io.outVerbose << bestVar << ",";
    *io.outVerbose << depth;
    *io.outVerbose << std::endl;
    //varImp[bestVar] += treeImprovement;
#endif

    if (par.yaimp_flag) {
      //		  yaimp->add(bestVar, bestVar, treeImprovement);
      //		  if ((parent != NULL) && (parent->varID != bestVar))
      //        yaimp->add(parent->varID, bestVar, treeImprovement + lastImprovement);
    }

    Helper::getClassifiedVectors<T>(data, rowMaskVec, classifier, classVec,
        classMaskMultiVec);

    // build node
    iNode = new INode<T> (classifier);
    iNode->resize(classifier->resultingNodeSize());
    iNode->setParent(parent);

		TIMEPROF_STOP("JTreeCtrl~~makeNodeIntern");

    for (i = 0; i < classMaskMultiVec.size(); ++i) {
#ifdef JTREEDEBUG
      *io.outVerbose << "classMaskMultiVec: " << classMaskMultiVec.size();
      *io.outVerbose << "classMaskMultiVec[" << i << "]";
      Helper::printVec(*(classMaskMultiVec[i]));
      *io.outVerbose << "classifier: " << *classifier << std::endl;
      *io.outVerbose << std::endl << std::endl;
#endif
      iNode->kids[classVec[i]] = makeNode(data, par, *(classMaskMultiVec[i]),
          trainSetWeight, depth + 1, iNode, treeImprovement);
      delete classMaskMultiVec[i];
    }

    // fill NULL-nodes with the most common value
    for (i = 0; i < iNode->kids.size(); ++i) {
      if (iNode->kids[i] == NULL) {
        iNode->kids[i] = new INode<T>(
						gen.fct.getTermResult(data, par, par.depVar, rowMaskVec, 
							colMaskVec, trainSetWeight), par.ncol);
      }
    }

    return iNode;
  }

  virtual Tree<T, uli_t> *makeTree(
      DataFrame<T> *data,
      std::vector<uli_t> *colMaskVec, //if NULL then take all variables
      uli_t treeId,
      bool downsampling_flag,
      //OUTPUT
      //       var.  , vector with tree IDs
      std::vector<uli_t> &trainDataRows, DataTreeSet &oobSet,
      std::vector<uli_t> &oobDataRows, DataFrame<double> *yaimp_ = NULL) {

    try {
      if (data == NULL)
        throw Exception(ERRORCODE_5);
      this->setData(data, colMaskVec);
      if (data->getncol() <= par.depVar)
        throw Exception(ERRORCODE_18);
      this->sizeCutoff = par.targetPartitionSize;
      uli_t i, row;
      INode<T> *treeRoot;
      std::vector<uli_t> maskVec;
      std::vector<uli_t> rndSample;

      this->trainSetWeight.clear();
      this->trainSetWeight.resize(getData()->getnrow(), 0.0);

      if (importance != NULL)
        delete importance;
      this->importance = gen.fct.newImportanceObject();
      this->importance->setBaseData(&par, &io, data);
      this->importance->init();
      this->importance->reset();

      // create yaimp matrix
      //      if (par.yaimp_flag) yaimp = yaimp_;

      // create varImp
      varImp.assign(par.ncol, 0);

			// set class weights
			std::vector<double > classwt(data->getnrow(), 1.0);
			const std::vector<uli_t > &varCat = data->varCategories[par.depVar];

			if (strcmp(par.classweights, "") != 0)
			{
				std::vector<double > classwtPerIdx;
				// split input string
				std::string strClasswt = std::string(par.classweights);
				size_t posBegin = 0, posEnd = strClasswt.find(";");
				while (posEnd != std::string::npos)
				{
					double wt;
				  std::stringstream(strClasswt.substr(posBegin, posEnd - posBegin)) >> wt;
				  classwtPerIdx.push_back(wt);
					posBegin = posEnd + 1;
					posEnd = strClasswt.find(";", posEnd + 1);

				  //*io.outVerbose << "Setting weight of class "
					//	<< data->varClassIndexer[par.depVar][classwtPerIdx.size()-1]
					//	<< " to "
					//	<< wt << std::endl;
				}

      	for (size_t i = 0; i < classwt.size(); ++i) {
					const uli_t &cl = data->depIdxVec[i];
					classwt[i] = classwtPerIdx[cl];
				}
			}
			
			//if (par.autoclassreweight_flag) {
			//	double nMax = (double)(*std::max_element(varCat.begin(), varCat.end()));
			//	//std::cout << nMax << std::endl;
      //	for (i = 0; i < classwt.size(); ++i) {
			//		const uli_t &cl = data->depIdxVec[i];
			//		classwt[i] = nMax / (double)varCat[(size_t) cl]; 
			//	  //std::cout << classwt[i] << std::endl;
			//	}
			//}

      // get sample weights (bagging)
      double popratio = 1;
			i = 0; 
			while (i < this->trainSetSize) {
        row = gsl_rng_uniform_int(par.rng, getData()->getnrow());
        //std::cout << "row " << row << std::endl;
        //row = (uli_t)floor(gsl_rng_uniform(rng) * getData()->getnrow());

        popratio = 1.0;

        // simulate equal class weighting
        if (par.weightsim_flag) {
				  const uli_t &cl = data->depIdxVec[row];
          popratio -= varCat[(size_t) cl] / (double) data->getnrow();
        }

        if (popratio > gsl_rng_uniform(par.rng)) {
          // no replacement?
          if (downsampling_flag)
            this->trainSetWeight[row] = classwt[row];
          else
            this->trainSetWeight[row] += classwt[row];

					++i;
        }
      }
      trainDataRows.clear();
      oobDataRows.clear();
			//Helper::printVec(this->trainSetWeight);

      // create training set / oob set
      for (i = 0; i < getData()->getnrow(); ++i) {

        if (this->trainSetWeight[i] != 0) {
          trainDataRows.push_back(i);
        } else {
          oobDataRows.push_back(i);
          oobSet.at(i, treeId, true);
        }
      }

      for (i = 0; i < trainDataRows.size(); ++i) {
        //std::cout << trainDataRows[i] << ",";
      }

      fitFctPar = new FittingFctPar<T> (NULL, data, par.depVar, &trainDataRows,
          colMaskVec, &trainSetWeight, 0, 0, 0, par.rng,
          //OUTPUT
          importance, NULL, NULL, false, &twoWayInteraction, NULL);

      std::vector<T> outVec;
      getData()->getRowMaskedCol(par.depVar, trainDataRows, outVec);

      treeRoot
          = makeNode(*data, par, trainDataRows, trainSetWeight, 0, NULL, 0);

      delete fitFctPar;

      return (new Tree<T, uli_t> (treeRoot));

    } catch (std::exception &e) {
      std::string *out = new std::string("JTreeCtrl::makeTree:");
      out->append(e.what());
      throw Exception(out->c_str());
    }
  }
  ;

  void getFreqFormNode( //INPUT
      INode<T> *inode,
      //OUTPUT
      std::map<uli_t, double> *freqMap) {

    //		(*freqMap)[inode->getVarID()]++;
    //
    //		for (uli_t i = 0; i < inode->size(); ++i) {
    //			if (inode->at(i) != NULL)
    //				this->getFreqFormNode(	static_cast<INode<T > *>(inode->at(i)),
    //										freqMap);
    //		}
  }

  void getFreqImportance(Tree<T, uli_t> *tree, DataFrame<T> *data, std::map<
      uli_t, double> *freqMap) {

    if (data == NULL)
      throw Exception(ERRORCODE_5);
    if (tree == NULL)
      throw Exception(ERRORCODE_10);
    if (tree->getRoot() == NULL)
      throw Exception(ERRORCODE_11);

    this->getFreqFormNode(static_cast<INode<T> *> (tree->getRoot()), freqMap);

  }

  inline Importance<T>* getVarImp() {
    return importance;
  }
  ;

private:
  void setData(DataFrame<T> *data, std::vector<uli_t> *idxVec = NULL) {

    this->data = data;
    this->trainSetSize = (uli_t) (getData()->getnrow());
    this->makeBuffer(idxVec);
  }
  ;
  inline DataFrame<T> *getData() const {
    return this->data;
  }
  ;
  DataFrame<T> *data;

  void makeBuffer(std::vector<uli_t> *idxVec) {
    uli_t i;

    if (idxVec == NULL) {
      if (colBuffer != NULL)
        delete[] colBuffer;
      if (rowBuffer != NULL)
        delete[] rowBuffer;

      colBuffer = new uli_t[getData()->getncol() - 1];
      rowBuffer = new uli_t[getData()->getnrow()];

      for (i = 0; i < getData()->getncol(); ++i)
        if (i != par.depVar)
          colBuffer[i - ((i < par.depVar) ? 0 : 1)] = i;

      for (i = 0; i < getData()->getnrow(); ++i)
        rowBuffer[i] = i;

      this->colBufferSize = getData()->getncol() - 1;
    } else {
      if (colBuffer != NULL)
        delete[] colBuffer;
      if (rowBuffer != NULL)
        delete[] rowBuffer;

      colBuffer = new uli_t[idxVec->size()];
      rowBuffer = new uli_t[getData()->getnrow()];

      for (i = 0; i < idxVec->size(); ++i)
        colBuffer[i] = idxVec->at(i);

      for (i = 0; i < getData()->getnrow(); ++i)
        rowBuffer[i] = i;

      this->colBufferSize = idxVec->size();
    }

    makeChoose();
  }
  inline uli_t *getColBuffer() const {
    return this->colBuffer;
  }
  ;
  inline uli_t *getRowBuffer() const {
    return this->rowBuffer;
  }
  ;
  uli_t *colBuffer;
  uli_t *rowBuffer;

  void makeChoose() {
    if (colChoose != NULL)
      delete[] colChoose;
    if (rowChoose != NULL)
      delete[] rowChoose;

    colChoose = new uli_t[par.mtry];
    rowChoose = new uli_t[this->trainSetSize];
  }
  inline uli_t *getColChooseBuffer() const {
    return this->colChoose;
  }
  ;
  inline uli_t *getRowChooseBuffer() const {
    return this->rowChoose;
  }
  ;
  uli_t *colChoose;
  uli_t *rowChoose;

  inline void setFreqMap(std::map<uli_t, T> *freqMap) {
    this->freqMap = freqMap;
  }
  inline std::map<uli_t, T> *getFreqMap() const {
    return this->freqMap;
  }
  ;
  std::map<uli_t, T> *freqMap;
  //
  // 	inline void setGetBestVariableFct(GETBESTVARFCT_ASVAR(T)) {
  // 		this->getBestVariable = getBestVariable;
  // 	}
  // 	GETBESTVARFCT_ASVAR(T);
  //
  // 	inline void setGetNodePerformanceFct(GETNODEPERFFCT_ASVAR(T)) {
  // 		this->getNodePerformance = getNodePerformance;
  // 	}
  // 	GETNODEPERFFCT_ASVAR(T);
  //
  // 	inline void setGetTermResultFct(GETTERMRESFCT_ASVAR(T)) {
  // 		this->getTermResult = getTermResult;
  // 	}
  // 	GETTERMRESFCT_ASVAR(T);
  //
  std::vector<double> trainSetWeight;
  uli_t trainSetSize;
  uli_t sizeCutoff;
  uli_t colBufferSize;

  Importance<T> *importance;
  std::vector<double> varImp;
  std::vector<std::pair<double, std::pair<uli_t, uli_t> > > twoWayInteraction;
  //	DataFrame<double > *yaimp;

  FittingFctPar<T> *fitFctPar;

  RJunglePar par;
  RJungleIO io;
  RJungleGen<T> gen;
};

#endif /*JTREECTRL_H_*/
