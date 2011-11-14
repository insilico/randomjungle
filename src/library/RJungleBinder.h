/*
 * RJungleBinder.h
 *
 *  Created on: 13.11.2008
 *      Author: schwarz
 */

#ifndef RJUNGLEBINDER_H_
#define RJUNGLEBINDER_H_

#include "Prediction.h"
#include "Importance.h"
#include "TImportance.h"
#include "PermImportance.h"
#include "Helper.h"
#include "CmplFct.h"
#include "Proximities.h"

template<class T>
class RJungleBinder {
public:
  RJungleBinder() :
    intrinsicImportance(NULL), permImportance(NULL), yaimp(NULL) {
  }
  ;

  virtual ~RJungleBinder() {
    if (trees.size() != 0) {
      typename std::vector<Tree<T, uli_t> *>::iterator it;
      it = trees.begin();
      while (it != trees.end()) {
        if (*it != NULL)
          delete *it;
        ++it;
      }
      trees.clear();
    }

    if (cmpldTrees.size() != 0) {
      typename std::vector<CmpldTree<T> *>::iterator it;
      it = cmpldTrees.begin();
      while (it != cmpldTrees.end()) {
        if (*it != NULL)
          delete *it;
        ++it;
      }
      cmpldTrees.clear();
    }

    if (this->yaimp != NULL)
      delete this->yaimp;

    if (this->intrinsicImportance != NULL)
      delete intrinsicImportance;

    if (this->permImportance != NULL)
      delete permImportance;
  }

  RJunglePar par;
  RJungleIO io;
  RJungleGen<T> gen;

  std::vector<std::pair<double, uli_t> > freqPairs;
  std::vector<std::pair<double, std::pair<uli_t, uli_t> > > twoWayInteraction;

  // cache variables
  std::vector<Tree<T, uli_t> *> trees;
  std::vector<CmpldTree<T> *> cmpldTrees;
  Importance<T> *intrinsicImportance;
  PermImportance<T> *permImportance;
  Prediction<T> pred;

  // out of bag data set
  DataTreeSet oobSet;

  // yast another importance
  DataFrame<double> *yaimp;

  // sample proximities
  Proximities<T> proxi;

  // variable proximity
  Proximities<T> varProxi;
};

#endif /* RJUNGLEBINDER_H_ */
