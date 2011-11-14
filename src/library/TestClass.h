#ifndef TESTCLASS_H_
#define TESTCLASS_H_


#include <iostream>
#include <sys/time.h>
#include <limits>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <cassert>
#include <sstream>

#include "treedefs.h"
#include "DataFrame.h"
#include "Helper.h"
#include "FittingFct.h"
#include "TMClassAtom.h"
#include "SClassAtom.h"
#include "Tree.h"
#include "JTreeCtrl.h"
#include "RJungleCtrl.h"

#include "lr.h"

class TestClass
{
public:
	TestClass();
	virtual ~TestClass();
	void testT2ClassAtom(void) const;
	void testTMClassAtom(void) const;
  void testSClassAtom(void) const;
	void testNode(void) const;
	void testDataFrame(void) const;
	void testDataFrameAndNode(void) const;
	void testFind(void) const;
	void testPurity() const;
	void testGini1() const;
	void testGini2() const;
	void testTree() const;
	void testJTreeCtrl4() const;
	void testJTreeCtrl5() const;
	void testRandomJungleCtrl1() const;
	void testRandomJungleCtrl2() const;
	void testCI() const;
	void testMvt() const;
	void testC_maxabsConditionalPvalue() const;
	void helpers(void) const;
	void lr(void) const;

	void train_save( char *savename, lr_predict *lrp) const;
	lr_predict *mk_train_lr_predict( spardat *sp, dym *factors, dyv *outputs, lr_options *opts) const;
};

#endif /*TESTCLASS_H_*/
