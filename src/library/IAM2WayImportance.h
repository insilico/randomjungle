/*
 * IAM2WayImportance.h
 *
 *  Created on: 02.02.2009
 *      Author: banoslo
 */

#ifndef IAM2WAYIMPORTANCE_H_
#define IAM2WAYIMPORTANCE_H_

#include "treedefs.h"
#include "Helper.h"
#include "Importance.h"

template <class T>
class IAM2WayImportance: public Importance<T> {
public:
  IAM2WayImportance(uli_t vecSizeLimit = 1000000) {
    this->vecSizeLimit = vecSizeLimit;
  }
  virtual ~IAM2WayImportance() {
  }

  virtual void print() {
    // create pairs
    Helper::vecToPairs(freqPairs, amounts);

    // sort intrinsicVar
    sort(freqPairs.begin(), freqPairs.end());

    //Helper::printVec(amounts);

    //for (uli_t i = 0; i < indexer.size(); ++i) {
    //  std::cout << indexer[i].first << " " << indexer[i].second << std::endl;
    //}

    Helper::printTreePair2Way<T>(
      freqPairs, indexer, *this->data, this->par->ntree, this->iteration,
      *this->io->outImportance);
  }

  inline void add(uli_t varID1, uli_t varID2, double val) {
    std::pair<uli_t, uli_t> idxPair
      = std::make_pair<uli_t, uli_t>(std::min(varID1, varID2), std::max(varID1, varID2));
    // find index
    idxIt = find(indexer.begin(), indexer.end(), idxPair);

    // indexes found?
    if (idxIt == indexer.end()) {
      // if pair was not found, check: index vector overflow?
      if (indexer.size() == vecSizeLimit) {
        // if overflow, delete pair with lowest score or do not insert new pair
        Helper::vecToPairs(freqPairs, amounts);

        sort(freqPairs.begin(), freqPairs.end());

        if (freqPairs.begin()->first < val) {
          idx = freqPairs.begin()->second;
          amounts[idx] = val;
          indexer[idx].first = idxPair.first;
          indexer[idx].second = idxPair.second;
        }
      } else {
        // if no overflow, push content back
        indexer.push_back(idxPair);
        amounts.push_back(val);
      }
    } else {
      // if pair was found, add value to it
      amounts[idxIt - indexer.begin()] += val;
    }
  }

  virtual void add(Importance<T> *impx) {
    IAM2WayImportance<T> *imp = (IAM2WayImportance<T>*) impx;

    for (uli_t i = 0; i < imp->indexer.size(); ++i) {
      add(imp->indexer[i].first, imp->indexer[i].second, imp->amounts[i]);
    }
  }

  virtual void reset() {
    amounts.clear();
  }

  static Importance<T>* newImportanceObject() {
    return new IAM2WayImportance<T> ();
  }


  uli_t vecSizeLimit;
  uli_t idx;
  std::vector<std::pair<uli_t, uli_t> >::iterator idxIt;
  std::vector<std::pair<double, uli_t> > freqPairs;
  std::vector<std::pair<uli_t, uli_t> > indexer;
  std::vector<double> amounts;
};

#endif /* IAM2WAYIMPORTANCE_H_ */
