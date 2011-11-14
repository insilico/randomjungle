/*
 * Importance.h
 *
 *  Created on: 29.01.2009
 *      Author: schwarz
 */

#ifndef IMPORTANCE_H_
#define IMPORTANCE_H_

#include <vector>

#include "treedefs.h"
#include "RJunglePar.h"
#include "RJungleIO.h"


template <class T >
class Importance {
public:
  Importance() {};

  virtual ~Importance() {};

  virtual void printHeader() {}
  virtual void print() {}
  virtual void save() {}

  void setBaseData(RJunglePar *par, RJungleIO *io, DataFrame<T > *data) {
    this->par = par;
    this->io = io;
    this->data = data;
  }

  void setBaseDataFrom(Importance<T > &imp) {
    this->par = imp.par;
    this->io = imp.io;
    this->data = imp.data;
  }

  virtual void add(Importance<T > *imp) {}

  virtual void init() {}
  virtual void reset() {}

  virtual void setIteration(uli_t iteration) {this->iteration = iteration;};
  virtual uli_t getIteration() {return iteration;};

  virtual void getBestVariables(size_t size, std::vector<uli_t > &outVec) {};
  virtual void combineMpi() {};

  uli_t iteration;
  RJunglePar *par;
  RJungleIO *io;
  DataFrame<T > *data;

};

#endif /* IMPORTANCE_H_ */
