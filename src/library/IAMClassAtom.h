/*
 * IAMClassAtom.h
 *
 *  Created on: 27.01.2009
 *      Author: schwarz
 */

#ifndef IAMCLASSATOM_H_
#define IAMCLASSATOM_H_

#include <boost/dynamic_bitset.hpp>
#include "ClassAtom.h"
#include "Helper.h"

#ifndef Tvector
#define Tvector std::vector<T >
#endif
#ifndef CCvector
#define CCvector std::vector<uli_t >
#endif

template <class T>
class IAMClassAtom: public ClassAtom<T, uli_t > {
public:
  IAMClassAtom() { this->setErrorMissing(2);};

  virtual ~IAMClassAtom(void) { };

  /*
   * Classifies an input
   */
  // TODO: function is wrong
  virtual uli_t classify(Tvector &sample) const {
    uli_t i;

    for (i = 0; i < varID.size(); ++i) {
      if (sample[varID[i]] == this->getMissingCode())
        return this->getErrorMissing();
    }

    uli_t counter;
    boost::dynamic_bitset<> isLeft;
    isLeft.resize(varID.size(), false);

    for (i = 0; i < varID.size(); ++i)
      for (counter = 0; counter < classSet[i].size(); ++counter)
        if (sample[varID[i]] == classSet[i][counter]) {
          isLeft.set(i, true);
          break;
        }

    return (isLeft.any())?0:1;
  }

  // TODO: function is wrong
  virtual uli_t classify(T *sample) const {
    uli_t i;

    for (i = 0; i < varID.size(); ++i) {
      if (rjungle::magicAt(sample, varID[i]) == this->getMissingCode())
        return this->getErrorMissing();
    }

    uli_t counter;
    boost::dynamic_bitset<> isLeft;
    isLeft.resize(varID.size(), false);

    for (i = 0; i < varID.size(); ++i)
      for (counter = 0; counter < classSet[i].size(); ++counter)
        if (rjungle::magicAt(sample, varID[i]) == classSet[i][counter])
          isLeft.set(i, true);

    return (isLeft.any())?0:1;
  }

  bool isConsistent() const {
    return true;
  }

  virtual std::ostream& print(std::ostream& os) const {
    uli_t i;
    os << "iamclass:";
    os << "varID:";
    Helper::printVec<uli_t >(varID);
    os << "classSet:";
    for (i = 0; i < classSet.size(); ++i) {
      os << ",";
      Helper::printVec<T >(classSet[i]);
    }
    os << std::endl;
    return os;
  }

  virtual void printXml(std::ostream& os) const {
    /*
    os
    << "<classifier id=\"" << ca_IAM
    << "\" size=\"" << 0
    << "\" varID=\"" << varID
    << "\" type=\"sclass\">Not impl. yet</classifier>" << std::endl;
    */
  }

  inline const std::vector<std::vector<T > > &getClassSet() const {return classSet;}
  inline void setClassSet(std::vector<std::vector<T > > &set) {this->classSet = set;}

  virtual inline uli_t resultingNodeSize() const { return 2; }

  inline bool operator==(const IAMClassAtom<T> &iAMClassAtom) {
    return iAMClassAtom.classSet == this->classSet;
  }

  virtual inline uli_t getType() const { return ca_IAM; }

  inline void setVarID(std::vector<uli_t > &varID) {this->varID = varID; }
  inline std::vector<uli_t > getVarID() const {return this->varID;}

  virtual bool isThereVarID(uli_t varID) const {
    for (uli_t i = 0; i < this->varID.size(); ++i)
      if (varID == this->varID[i]) return true;

    return false;
  }

  std::vector<uli_t > varID;
  std::vector<std::vector<T > > classSet;
};

/*
 * Source
 */

// ostream
template <class T>
std::ostream& operator<<(std::ostream& os, const IAMClassAtom<T> &iAMClassAtom) {
  return iAMClassAtom.print(os);
}


#endif /* IAMCLASSATOM_H_ */
