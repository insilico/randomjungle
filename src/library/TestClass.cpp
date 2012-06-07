/*
 * Includes
 */

#include <vector>
#include <string>

#include "TestClass.h"
#include "CondInf.h"
#include "mvt.h"
#include "CondInf.h"
#include "SaveCollector.h"
#include "PermImportance.h"
#include "Profiler.h"
#include "TimeProf.h"
#include "librjungle.h"

#include "lr.h"

/*
 * Source
 */
TestClass::TestClass() { }

TestClass::~TestClass() {
}

/** 
 * \brief Test logistic regression
 */
void TestClass::train_save( char *savename, lr_predict *lrp) const {
	PFILE *f;
	f = sure_pfopen( savename, "w");
	out_lr_predict( f, lrp);
	sure_pfclose( f, savename);
}

lr_predict *TestClass::mk_train_lr_predict( spardat *sp, dym *factors, dyv *outputs,
		lr_options *opts) const{

	lr_train *lrt;
	lr_predict *lrp;
	lrt = mk_lr_train( sp, factors, outputs, NULL, opts);
	lrp = mk_lr_predict( lrt_b0_ref(lrt), lrt_b_ref(lrt));
	free_lr_train( lrt);
	return lrp;
}

void TestClass::lr(void) const {

  time_t start, stop;
  dyv *outputs;
  dym *factors;
  spardat *sp;
  lr_options *opts;
  lr_predict *lrp;

  /* Parse lr options. */
  opts = mk_lr_options();

  /* Load data. */
  sp = NULL;
  factors = NULL;
  outputs = NULL;

	mk_read_dym_for_csv((char *) "./tests/liblrtest.csv", &factors, &outputs);
	if (!dyv_is_binary( outputs)) {
		my_error("run_train: Error: csv output column is not binary.\n");
	}

  /* Train. */
  start = time( NULL);
  lrp = mk_train_lr_predict( sp, factors, outputs, opts);
  stop = time( NULL);

  /* Clean. */
  if (sp != NULL) free_spardat( sp);
  if (factors != NULL) free_dym( factors);
  if (outputs != NULL) free_dyv( outputs);
  free_lr_options( opts);

  /* Save enough info to make predictions later. */
  train_save((char *) "./tests/liblrtest.output", lrp);
  free_lr_predict( lrp);

  /* Output stuff if desired. */
  if (Verbosity >= 1) {
    printf( "TOTAL ALG TIME = %f seconds\n", difftime(stop, start));
  }
  
}

void TestClass::helpers(void) const {

  // random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, 12345);

  // declare variables
  std::vector<double> vec, vec2;

  // fill vector
  vec.clear();
  vec.push_back(-1.1);
  vec.push_back(2.2);
  vec.push_back(33.33);

  // median
  assert(Helper::getMedian<double>(vec) == 2.2);

  // diff
  vec2 = Helper::getDiff<double>(vec, 1.1);
  assert(vec2[0] == -2.2);
  assert(vec2[1] == 1.1);
  assert(vec2[2] == 32.23);

  // median
  vec.push_back(444.444);
  assert(ceil(Helper::getMedian<double>(vec)) == 18);

  // mean
  assert(floor(Helper::getMean<double>(vec) * 1e+04) == 1197185);

  // abs
  assert((Helper::getAbs<double>(vec)).at(0) == 1.1);

  // mad
  assert(floor(Helper::getMad<double>(vec) * 1e+06) == 25522959);

  // get ranks
  std::vector<size_t> vec3;
  vec3 = Helper::getRanks<double>(vec2);
  assert(vec3[0] == 0);
  assert(vec3[1] == 1);
  assert(vec3[2] == 2);

  vec3 = Helper::getRanks<double>(vec2, false);
  assert(vec3[0] == 2);
  assert(vec3[1] == 1);
  assert(vec3[2] == 0);

  // save data collector
  SaveCollector saveCol;

  std::vector<double> doubleVec1, doubleVec2;
  std::vector<int> intVec1;
  std::vector<size_t> size_tVec1;
  std::vector<uli_t> uli_tVec1;
  std::vector<std::string> stringVec1;

  doubleVec1.push_back(1.1);
  doubleVec1.push_back(1.2);
  doubleVec1.push_back(1.3);

  intVec1.push_back(1);
  intVec1.push_back(2);
  intVec1.push_back(3);

  doubleVec2.push_back(3.1);
  doubleVec2.push_back(3.2);
  doubleVec2.push_back(3.3);

  size_tVec1.push_back(10);

  uli_tVec1.push_back(100);
  uli_tVec1.push_back(200);

  stringVec1.push_back("five");
  stringVec1.push_back("six");
  stringVec1.push_back("eins");

  saveCol.push_back(&doubleVec1, "Double1");
  saveCol.push_back(&intVec1, "Int1");
  saveCol.push_back(&doubleVec2, "DoubleTwo");
  //saveCol.push_back(&size_tVec1, "size_t");
  //saveCol.push_back(&uli_tVec1, "uli_t");
  saveCol.push_back(&stringVec1, "A_String");

  //std::cout << saveCol << std::endl;

  saveCol.orderRow.push_back(0);
  saveCol.orderRow.push_back(1);
  saveCol.orderRow.push_back(2);

  //std::cout << saveCol << std::endl;

  saveCol.orderRow.clear();
  saveCol.orderRow.push_back(2);
  saveCol.orderRow.push_back(0);
  saveCol.orderRow.push_back(1);

  //std::cout << saveCol << std::endl;

  // permutation importance
  // init objects
#ifndef __SPARSE_DATA__
  RJungleGen<double> gen;
  PermImportance<double> permImp;
  DataFrame<double> *data;
  RJunglePar par;
  RJungleIO io;
  // set parameters
  par = initRJunglePar();
  par.filename = (char *) "./tests/data_from_r.dat";
  par.outprefix = (char *) "./tests/testlib_helpers";
  par.depVarName = (char *) "Species";
  par.rng = r;
  par.ntree = 1;
  par.mtry = 2;
  par.seed = 123;
  io.open(par);

  // get data
  Helper::getCSVdim(par);
  data = new DataFrame<double> (par);
  par = data->getDataFromCSV(); // and return updated parameters
  gen.init(par, *data);

  // check seq
  std::vector<uli_t> idx, idx2;
  idx2.push_back(0);
  idx2.push_back(1);
  Helper::seq<uli_t>(idx, 0, 1, 1);
  assert(idx == idx2);

  // test pearson's correlation coefficient
  Helper::seq<uli_t>(idx, 0, par.nrow - 1, 1);
  assert(floor(Helper::getPearsonCor<double>(*data, idx, 0, 1) * 1e+07)
      == -1175698);

  // Init tree controller
  JTreeCtrl<double> treeCtrl;
  treeCtrl = JTreeCtrl<double> (par, io, gen);

  // Build tree initialization
  std::vector<uli_t> trainDataRows;
  std::vector<uli_t> oobDataRows;
  DataTreeSet oobSet;
  DataFrame<double> *yaimp = NULL;
  Tree<double, uli_t> *onetree = NULL;
  CmpldTree<double> *oneCmpldtree = NULL;
  std::vector<uli_t> *colMaskVec = NULL;

  trainDataRows.clear();
  oobDataRows.clear();
  oobSet.init(par.nrow, par.ntree);

  // build tree
  onetree = treeCtrl.makeTree(data, colMaskVec, 0, par.downsampling_flag,
      trainDataRows, oobSet, oobDataRows, yaimp);
  //std::cout << *onetree << std::endl;
  oneCmpldtree = gen.fct.cmplTree(*onetree);

  // get cutoffs
  std::vector<double> cutoffs, cutoffs2;
  gen.fct.getCutoffs(*oneCmpldtree, 1, cutoffs);
  gen.fct.getCutoffs(*oneCmpldtree, 2, cutoffs2);
  //Helper::printVec(cutoffs);
  //Helper::printVec(cutoffs2);

  // Interval assignment
  std::vector<uli_t> groups, groups2;

  cutoffs.clear();
  cutoffs.push_back(2.45);
  cutoffs.push_back(2.55);

  cutoffs2.clear();
  cutoffs2.push_back(2.45);
  cutoffs2.push_back(4.7);
  cutoffs2.push_back(4.85);
  cutoffs2.push_back(5.05);

  Helper::seq<uli_t>(idx, 0, par.nrow - 1, 1);
  Helper::getIntervalGroup(*data, idx, 1, cutoffs, groups);
  Helper::getIntervalGroup(*data, idx, 2, cutoffs2, groups2);

  assert(groups[0] == 2);
  assert(groups[100] == 2);
  assert(groups[41] == 0);
  assert(groups[53] == 0);
  assert(groups[119] == 0);

  assert(groups2[0] == 0);
  assert(groups2[100] == 4);
  assert(groups2[41] == 0);
  assert(groups2[53] == 1);
  assert(groups2[119] == 3);

  // create conditional importance grid
  std::vector<uli_t> grid;
  Helper::addToCIGrid(groups, grid);
  Helper::addToCIGrid(groups2, grid);

  assert(grid[0] == 1);
  assert(grid[100] == 10);
  assert(grid[41] == 0);
  assert(grid[53] == 2);
  assert(grid[119] == 6);

  // permute in groups
  Helper::seq<uli_t>(idx, 0, par.nrow - 1, 20);
  Helper::getIntervalGroup(*data, idx, 1, cutoffs, groups);
  Helper::getIntervalGroup(*data, idx, 2, cutoffs2, groups2);
  grid.clear();
  Helper::addToCIGrid(groups, grid);
  Helper::addToCIGrid(groups2, grid);

  double *vecPtr = new double[idx.size()];

  data->getColMaskedRow(1, idx, vecPtr);
  Helper::permuteInGrid<double>(r, vecPtr, grid);

  assert(vecPtr[0] == 3.4);
  assert(vecPtr[1] == 3.5);
  assert(vecPtr[2] == 3.5);
  assert(vecPtr[3] == 2.4);
  assert(vecPtr[4] == 2);
  assert(vecPtr[5] == 3.3);
  assert(vecPtr[6] == 3.2);
  assert(vecPtr[7] == 3.1);

  delete[] vecPtr;
  // setBaseData
  //permImp.setBaseData(&par, &io, data);
  /*
   permImp.init();

   // set random data
   permImp.add(0, 1.2);
   permImp.add(1, 11.1);
   permImp.add(2, 0.1);
   permImp.add(3, -0.1);

   // output
   permImp.enableHeader();

   permImp.setIteration(0);
   permImp.save();

   // set random data
   permImp.add(1, 31.1);
   permImp.add(2, 122.1);
   permImp.add(3, -12.1);

   // output
   permImp.disableHeader();
   permImp.setIteration(1);
   permImp.save();
   */
	
	// test the profiler
	Profiler pr1;
	pr1.start();
	pr1.stop();

	TimeProf tpr1;

	tpr1.start("prof1");
	tpr1.start("prof2");
	tpr1.stop("prof1");
	tpr1.stop("prof2");

	tpr1.start("prof1");
	tpr1.start("prof2");
	tpr1.stop("prof1");
	tpr1.stop("prof2");

  // clean up
  delete oneCmpldtree;
  delete onetree;
  delete data;
  io.close();
#endif
  // free gsl
  gsl_rng_free(r);
}

void TestClass::testC_maxabsConditionalPvalue(void) const {
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, 12345);

  // declare variables
  double ans;
	CondInf<double> ci;
  ci.rng = r;

  /* test C_maxabsConditionalPvalue function */
  /* test 1 */

  ci.dCov_T.clear();
  ci.dCov_T.push_back(1.0);
  ci.dCov_T.push_back(0.0);
  ci.dCov_T.push_back(0.0);
  ci.dCov_T.push_back(1.0);

  ci.cMax = 1.96;
  ci.p = 1;
  ci.q = 2;
  ans = ci.C_maxabsConditionalPvalue();
  assert(floor(ans * 1e+06) == 97492);

  /* test 2 */

  int i, j;
  ci.dCov_T.clear();
  ci.dCov_T.assign(16, 0.0);

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      ci.dCov_T[i * 4 + j] = (i == j) ? 1.0 : 0;

  ci.q = 4;
  ans = ci.C_maxabsConditionalPvalue();
  assert(floor(ans * 1e+06) == 185479);

  /* test 3 */
  ci.dCov_T.clear();

  ci.dCov_T.push_back(1.0);
  ci.dCov_T.push_back(0.1);
  ci.dCov_T.push_back(0.2);
  ci.dCov_T.push_back(0.1);
  ci.dCov_T.push_back(0.1);
  ci.dCov_T.push_back(1.0);
  ci.dCov_T.push_back(0.3);
  ci.dCov_T.push_back(0.2);
  ci.dCov_T.push_back(0.2);
  ci.dCov_T.push_back(0.3);
  ci.dCov_T.push_back(1.0);
  ci.dCov_T.push_back(0.3);
  ci.dCov_T.push_back(0.1);
  ci.dCov_T.push_back(0.2);
  ci.dCov_T.push_back(0.3);
  ci.dCov_T.push_back(1.0);

  ans = ci.C_maxabsConditionalPvalue();
  assert(floor(ans * 1e+03) == 179);

  /* free gsl */

  gsl_rng_free(r);
}

void TestClass::testMvt(void) const {
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, 12345);

  // declare variables
  int *infin;
  int n, nu, maxpts, inform;
  double *lower, *upper, *corr, *delta;
  double abseps, releps, error, value;

  // define test values
  n = 5;
  nu = 4;
  error = 0;
  value = 0;
  inform = 0;
  maxpts = 25000;
  abseps = 0.001;
  releps = 0.0;

  /*
   lower <- carg$lower/sqrt(diag(carg$sigma))
   upper <- carg$upper/sqrt(diag(carg$sigma))
   corr <- cov2cor(carg$sigma)
   */

  lower = new double[5];
  lower[0] = -1;
  lower[1] = -1;
  lower[2] = -1;
  lower[3] = -1;
  lower[4] = -1;

  upper = new double[5];
  upper[0] = 3;
  upper[1] = 3;
  upper[2] = 3;
  upper[3] = 3;
  upper[4] = 3;

  corr = new double[10];
  corr[0] = 0.5;
  corr[1] = 0.5;
  corr[2] = 0.5;
  corr[3] = 0.5;
  corr[4] = 0.5;
  corr[5] = 0.5;
  corr[6] = 0.5;
  corr[7] = 0.5;
  corr[8] = 0.5;
  corr[9] = 0.5;

  delta = new double[5];
  delta[0] = 0.0;
  delta[1] = 0.0;
  delta[2] = 0.0;
  delta[3] = 0.0;
  delta[4] = 0.0;

  infin = new int[5];
  infin[0] = 2;
  infin[1] = 2;
  infin[2] = 2;
  infin[3] = 2;
  infin[4] = 2;

  // call function
  mvtdst_(&n, &nu, lower, upper, infin, corr, delta, &maxpts, &abseps, &releps,
      &error, &value, &inform, &r);

  assert(floor(value * 1e+03) == 506);

  delete[] corr;
  delete[] lower;
  delete[] upper;
  delete[] infin;
  delete[] delta;
}

void TestClass::testCI(void) const {
  CondInf<char> ci;
  uli_t i = 0;
  double ans = 0;

  /*
   * Check initals
   */
  assert(ci.x.size() == 0);
  assert(ci.y.size() == 0);
  assert(ci.w.size() == 0);
  assert(ci.t.size() == 0);
  assert(ci.dExp_T.size() == 0);
  assert(ci.dCov_T.size() == 0);
  assert(ci.dExp_y.size() == 0);
  assert(ci.dCov_y.size() == 0);

  /*
   * Check sum of w and w * x
   */
  ci.p = 1; // x is numerical
  uli_t n = 10;

  ci.x.clear();
  for (i = 0; i < n; ++i)
    ci.x.push_back(i);

  ci.w = std::vector<double>(n, 1);
  ci.calcSumw();
  ci.calcSumwxCT1();

  assert(ci.sumw == n);
  assert(ci.swx[0] == 45);
  /*
   * Check y rank trafo. and set q
   */

  ci.y.clear();
  ci.y.push_back(3);
  ci.y.push_back(1);
  ci.y.push_back(10);
  ci.y.push_back(10);

  ci.trafoRank(ci.y);

  assert(ci.y[0] == 1);
  assert(ci.y[1] == 0);
  assert(ci.y[2] == 2);
  assert(ci.y[3] == 2);

  ci.setQ();
  assert(ci.q == 1);

  ci.yIsNominal = true;

  ci.setQ();
  assert(ci.q == 3);

  /*
   * Check expectation
   */

  // y = (0, 0, 0, 0, 0, 1, 1, 1, 1, 1)^T
  ci.y.clear();
  for (i = 0; i < n; ++i)
    ci.y.push_back((i < n / 2) ? 0 : 1);

  ci.calcdExp_y();

  // E = (0.5, 0.5)^T
  assert(ci.dExp_y.size() == ci.q);
  assert(ci.dExp_y[0] == 0.5);
  assert(ci.dExp_y[1] == 0.5);

  /*
   * Check cov.
   */
  n = 3;

  ci.y.clear();
  ci.y.push_back(0);
  ci.y.push_back(1);
  ci.y.push_back(2);

  ci.w = std::vector<double>(n, 1);

  ci.calcSumw();
  ci.trafoRank(ci.y);
  ci.yIsNominal = true;
  ci.setQ();
  ci.calcdExp_y();
  ci.calcdCov_y();

  // 2nd test
  n = 10;

  ci.y.clear();
  for (i = 0; i < n; ++i)
    ci.y.push_back((i < n / 2) ? 0 : 1);

  ci.w = std::vector<double>(n, 1);

  ci.calcSumw();
  ci.trafoRank(ci.y);
  ci.yIsNominal = true;
  ci.setQ();
  ci.calcdExp_y();
  ci.calcdCov_y();

  // V = ( 0.25, -0.25)
  //     (-0.25,  0.25)
  assert(ci.dCov_y.size() == ci.q * ci.q);
  assert(ci.dCov_y[0] == 0.25);
  assert(ci.dCov_y[1] == -0.25);
  assert(ci.dCov_y[2] == -0.25);
  assert(ci.dCov_y[3] == 0.25);

  /*
   * Check sigma
   */
  ci.x.clear();
  ci.x.push_back(0);
  ci.x.push_back(1);
  ci.x.push_back(2);
  ci.x.push_back(0);
  ci.x.push_back(1);
  ci.x.push_back(2);
  ci.x.push_back(0);
  ci.x.push_back(1);
  ci.x.push_back(2);
  ci.x.push_back(0);

  ci.trafoRank(ci.x);
  ci.xIsNominal = true;
  ci.setP();

  ci.trafoRank(ci.y);
  ci.yIsNominal = true;
  ci.setQ();

  ci.calcSumw();
  ci.calcSumwxCT1();
  ci.calcdExp_y();
  ci.calcdCov_y();
  ci.calcdCov_T();
  ci.calcdExp_T();
  ci.calcT();
  ci.calcCMax();

  /*
   # R-Code
   library(party)
   attach(asNamespace("party"))
   yf <- gl(2, 5)
   xf <- rep(gl(3, 1),4)[1:10]
   weights <- rep(1,10)
   x <- sapply(levels(xf), function(l) as.numeric(xf == l))
   colnames(x) <- NULL
   y <- sapply(levels(yf), function(l) as.numeric(yf == l))
   colnames(y) <- NULL
   LinearStatistic(x, y, weights)
   linstat <- LinearStatistic(x, y, weights)
   expcov <- ExpectCovarLinearStatistic(x, y, weights)
   maxabs <- max(abs(linstat - expcov@expectation) / sqrt(diag(expcov@covariance)))
   ceiling(maxabs * 100) == 66
   */

	// FIX: Sometimes this will fail (why?)
  //assert(ceil(ci.cMax * 100) == 66);

  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, 12345);
  ci.rng = r;

  ans = ci.C_maxabsConditionalPvalue();
  //assert(floor(ans * 1e+03) == 789);
	/*
	 * # R-Code (for x,y see above)
	 * res = .Call("R_IndependenceTest", x, y, weights,new("LinStatExpectCovar", ncol(x), ncol(y)), new("VariableControl"), PACKAGE = "party")
	 * tstat = res[1]
	 * pval = res[2]
	 * floor(pval * 1000) == 789
	 */

  gsl_rng_free(r);

  /*
   * test 2
   */

  ci.x.clear();
  ci.x.push_back(0);
  ci.x.push_back(0);
  ci.x.push_back(1);
  ci.x.push_back(1);
  ci.x.push_back(1);
  ci.x.push_back(1);
  ci.x.push_back(1);
  ci.x.push_back(1);
  ci.x.push_back(2);
  ci.x.push_back(2);

  ci.xIsNominal = false;
  ci.setP();

  ci.trafoRank(ci.y);
  ci.yIsNominal = true;
  ci.setQ();

  ci.calcSumw();
  ci.calcSumwxCT1();
  ci.calcdExp_y();
  ci.calcdCov_y();
  ci.calcdCov_T();

  /*
   * Check Mu
   */
  ci.calcdExp_T();

  assert(ci.dExp_T.size() == 2);
  assert(ci.dExp_T[0] == 5);
  assert(ci.dExp_T[1] == 5);

  /*
   * Check statistic
   */
  ci.calcT();

  assert(ci.t.size() == 2);
  assert(ci.t[0] == 3);
  assert(ci.t[1] == 7);

  /*
   * Calc. cMax
   */
  ci.calcCMax();

  //assert(round(ci.cMax * 10) == 19);

  /*
   * Calc. asym. p value
   */
  ci.calcPAsymCase1();

  //assert(round(ci.pAsym * 100) == 6);
	/* 
	 * # R-Code (for y see above)
	 * .Call(
	 *   "R_IndependenceTest", matrix(rank(c(0,0,1,1,1,1,1,1,2,2)), ncol = 1),
	 *   y, weights,new("LinStatExpectCovar", ncol(x), ncol(y)), new("VariableControl"),
	 *   PACKAGE = "party")
	 */

  /*
   * Calc. split
   */
  ci.x.clear();
  ci.x.push_back(3);
  ci.x.push_back(0);
  ci.x.push_back(0);
  ci.x.push_back(1);
  ci.x.push_back(2);
  ci.x.push_back(3);
  ci.x.push_back(3);
  ci.x.push_back(3);
  ci.x.push_back(3);
  ci.x.push_back(3);

  CondInf<char>::calcSplitCase1(ci);

  assert(ci.split == 2);

  /*
   * case 2 (two-sample)
   */
  ci.xIsNominal = false;
  ci.setP();
  ci.calcCMaxCase2();
  ci.calcPAsymCase2();

  //assert(floor(ci.cMax * 10) == 22);
  //assert(floor(ci.pAsym * 1000) == 26);

  /*
   * Calc. split case 2
   */
  CondInf<char>::calcSplitCase2(ci);

  assert(ci.split == 2);

}

void TestClass::testSClassAtom(void) const {
  /*
   // define testTClassAtom1 with a 0 sized vector
   SClassAtom<int > s_1;
   assert(s_1.size() == 0); //define with a 0 sized vector

   std::vector<int > testSet;
   testSet.push_back(1);
   s_1.setClassSet(testSet);
   assert(s_1.getClassSet() == testSet);

   assert(s_1.classify(0) == 1);
   assert(s_1.classify(1) == 0);
   assert(s_1.classify(2) == 1);

   // define with a 1 sized vector
   SClassAtom<int > s_2 = s_1;
   assert(s_2.getClassSet() == testSet);
   assert(s_2 == s_1);

   assert(s_2.resultingNodeSize() == 2);
   assert(s_2.getType() == ca_S);

   assert(s_2.classify(0) == 1);
   assert(s_2.classify(1) == 0);
   assert(s_2.classify(2) == 1);
   */
}

/*
 void TestClass::testT2ClassAtom(void) const
 {
 // define testTClassAtom1 with a 0 sized vector
 T2ClassAtom<int > t2_1;
 assert(t2_1.size() == 0); //define with a 0 sized vector
 t2_1.setThreshold(2);
 assert(t2_1.getThreshold() == 2);

 // define with a 1 sized vector
 T2ClassAtom<int > t2_2 = t2_1;
 assert(t2_2.getThreshold() == 2);
 assert(t2_2 == t2_1);

 assert(t2_2.resultingNodeSize() == 2);
 assert(t2_2.getType() == ca_T2);

 t2_2.setThreshold(1);
 assert(t2_2.getThreshold() == 1);
 assert(t2_2.classify(0) == 0);
 assert(t2_2.classify(1) == 0);
 assert(t2_2.classify(2) == 1);

 T2ClassAtom<char > t2_3;
 t2_3.setThreshold(2);
 assert(t2_3.getThreshold() == 2);

 T2ClassAtom<double > t2_4;
 t2_4.setThreshold(1.5);
 assert(t2_4.getThreshold() == 1.5);
 }

 void TestClass::testTMClassAtom(void) const {

 // define a 2 sized vector
 TMClassAtom<int > testTMClassAtom1;
 std::vector<int > intVec;
 intVec.push_back(1);
 intVec.push_back(5);
 intVec.push_back(2);
 testTMClassAtom1.setThresholds(intVec);

 std::vector<int > &intVec2 = testTMClassAtom1.getThresholds();
 assert(intVec2[0] == 1);
 assert(intVec2[1] == 2);
 assert(intVec2[2] == 5);

 assert(testTMClassAtom1.size() == 4);
 assert(testTMClassAtom1.resultingNodeSize() == 4);

 //checking stream operator
 std::ostringstream s1;
 s1 << testTMClassAtom1;
 std::string s2 = s1.str();
 assert(s2.compare("tmclass(<=):1,2,5") == 0);

 // checking test cases
 assert(testTMClassAtom1.classify(-1) == 0);
 assert(testTMClassAtom1.classify(0) == 0);
 assert(testTMClassAtom1.classify(1) == 0);
 assert(testTMClassAtom1.classify(2) == 1);
 assert(testTMClassAtom1.classify(3) == 2);
 assert(testTMClassAtom1.classify(4) == 2);
 assert(testTMClassAtom1.classify(5) == 2);
 assert(testTMClassAtom1.classify(6) == 3);
 }

 void TestClass::testNode(void) const {
 // no need to delete this ptr at end of function
 // this will be done by destructor of node2
 T2ClassAtom<double > *t2_1 = new T2ClassAtom<double >();
 t2_1->setThreshold(0.5);
 assert(t2_1->getThreshold() == 0.5);

 // using _late binding_ functionality! (upcast)
 // typeOf(*t2_1) = T2ClassAtom<double >
 // typeOf(*node2.classifier)   = T2ClassAtom<double > (base-class ClassAtom)
 INode<double > n1(t2_1, 11);
 assert(n1.varID == 11);

 // the same output, because c++ using magic "vptr"/"vtable"(late binding)
 std::ostringstream s1;
 s1 << *t2_1;
 std::string s1_ = s1.str();
 std::ostringstream s2;
 s2 << *(n1.getClassifier());
 std::string s2_ = s2.str();
 assert(s1_.compare(s2_) == 0);

 // test a small node split
 TermClassAtom<double, uli_t > *t2_1a;
 TermClassAtom<double, uli_t > *t2_1b;

 t2_1a = new TermClassAtom<double, uli_t >(-1.0);
 t2_1b = new TermClassAtom<double, uli_t >(1.0);

 n1.push_back(new INode<double >(t2_1a, 10));
 n1.push_back(new INode<double >(t2_1b, 20));

 assert(((INode<double > *)n1[0])->classifier == t2_1a);
 assert(((INode<double > *)n1[1])->classifier == t2_1b);

 //classify
 std::vector<double > smpl(30,0);
 assert(n1.classify(smpl) == -1.0);

 }

 void TestClass::testDataFrame(void) const
 {
 uli_t i, j;

 DataFrame<double > data(2, 3);

 for(i = 0; i < data.getnrow(); ++i)
 for(j = 0; j < data.getncol(); ++j)
 data.setArray(i, j, j + i * data.getnrow());
 for(i = 0; i < data.getnrow(); ++i)
 for(j = 0; j < data.getncol(); ++j)
 assert(data.getArray(i, j) == j + i * data.getnrow());

 assert(data.at(0, 0) == data.at(0, 0));

 std::vector<uli_t > maskVec;
 std::vector<double > outVec(3, 0);
 maskVec.push_back(1);
 maskVec.push_back(0);
 data.getRowMaskedCol(0, maskVec, outVec);
 assert(outVec[0] == 2);
 assert(outVec[1] == 0);
 maskVec.push_back(2);
 data.getColMaskedRow(0, maskVec, outVec);
 assert(outVec[0] == 1);
 assert(outVec[1] == 0);
 assert(outVec[2] == 2);

 }

 void TestClass::testFind(void) const { }

 void TestClass::testDataFrameAndNode(void) const { }
 */

void TestClass::testJTreeCtrl5() const {
  /*
   double missingcode = -999;
   uli_t ntree = 1;
   uli_t dataSet = 1;
   uli_t depVar = 0;
   uli_t i;
   uli_t nrow;
   uli_t ncol;

   if (dataSet == 1) {
   nrow = 2484;
   ncol = 7;
   } else if (dataSet == 2) {
   nrow = 3500;//3500;
   ncol = 9188;//9188;
   } else if (dataSet == 3) {
   //regression
   nrow = 10;
   ncol = 4;
   } else if (dataSet == 4) {
   depVar = 3622;//SNP6_153 0 cM from DR, increases RA risk only in women
   nrow = 3500;//3500;
   ncol = 9100;//9188;
   }  else if (dataSet == 5) {
   //nrow = 1200;
   nrow = 120; // with 120 real samples and } 1000 additional zerosamples
   //, so we avoid the strobel problem
   ncol = 6;    // funny , but true and i dont know why!
   } else if (dataSet == 6) {
   nrow = 10;
   ncol = 4;
   } else if (dataSet == 7) {
   nrow = 1500;//3500;
   ncol = 9000;//9198;
   } else if (dataSet == 8) {
   nrow = 2603;
   ncol = 26267;//26267 only good SNPs //45707 all SNPs;
   }

   DataFrame<double > *data = new DataFrame<double >(nrow, ncol);
   timeval start, end;
   std::vector<Tree<double, uli_t > * > trees;
   std::map<uli_t, double> *freqMap = new std::map<uli_t, double>();
   std::map<double, uli_t> mapInverse;
   std::map<uli_t, double> decGiniMap;
   std::vector<uli_t > keys;
   std::vector<uli_t >::iterator itKey;
   data->setMissingCode(missingcode);
   data->setDepVar(depVar);

   if (dataSet == 1) {
   std::pair<uli_t, uli_t > csvdim = Helper::getCSVdim("Ger5.csv");
   std::cout << csvdim.first << " " << csvdim.second << std::endl;
   data->getDataFromCSV(
   "Ger5.csv",
   nrow,
   ncol,
   0,
   0,
   0,
   0);
   //std::cout 	<< "data:\n" << *data << std::endl;
   } else if (dataSet == 2) {
   data->getDataFromCSV(
   "plink_chrALL_rep0001.csv",
   nrow,
   ncol,
   0,
   5,
   1,
   6//6
   );
   } else if (dataSet == 3) {
   data->getDataFromCSV(
   "testdatareg.csv",
   nrow,
   ncol,
   0,
   0,
   1,
   1
   );
   } else if (dataSet == 4) {
   data->getDataFromCSV(
   "plink_chrALL_rep0001.csv",
   nrow,
   ncol,
   0,
   5,
   1,
   4
   );
   } else if (dataSet == 5) {
   data->getDataFromCSV(
   "info.csv",
   nrow,
   ncol,
   0,
   0,
   0,
   0
   );
   } else if (dataSet == 6) {
   data->getDataFromCSV(
   "testdata.csv",
   nrow,
   ncol,
   0,
   0,
   1,
   1
   );
   } else if (dataSet == 7) {
   data->getDataFromCSV(
   "plinkpheno_chrALL_rep0001.csv",
   nrow,
   ncol,
   0,
   9,//9, = antiCCP (SNP18_269 controls effect of DR on anti-CCP and increases RA risk), 10 = lgM(QTL for SNP11_389),11 = Severty(QTL for SNP9_187/SNP9_192),
   2001,
   16
   );
   } else if (dataSet == 8) {
   data->getDataFromCSV(
   "illumina_50k.csv",
   nrow,
   ncol,
   0,
   5,
   1,
   6
   );
   }


   const std::vector<std::string > &varNames = data->getVarNames();
   std::cout << "dep. Var. name:" << varNames[depVar] << std::endl;

   if (dataSet == 0) {
   Helper::printVec<std::string >(varNames);
   std::cout << "categories of data:" << std::endl;
   std::map<double, uli_t >::iterator itCat;
   for (i = 0; i < data->getncol(); ++i) {
   std::cout << "var " << varNames[i] << ": ";
   itCat = data->varCategories[i].begin();
   while (itCat != data->varCategories[i].end()) {
   std::cout << itCat->first << " " << itCat->second << ",";
   ++itCat;
   }
   std::cout << std::endl;
   }
   }

   //std::cout << *data << std::endl;

   // define
   gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);

   JTreeCtrl<double > treeCtrl;

   GETBESTVARFCT_ASVAR(double);
   GETNODEPERFFCT_ASVAR(double);
   GETTERMRESFCT_ASVAR(double);

   if (dataSet == 3 || dataSet == 7) {
   getBestVariable = FittingFct::CARTreg<double >;
   getNodePerformance = Helper::puritySos<double >;
   getTermResult = Helper::getMean2<double >;
   } else {
   getBestVariable = FittingFct::CARTgini<double >;
   //		FittingFct::checkUsabilityLOTUSginiHybrid<double >(*data, depVar);
   //		getBestVariable = FittingFct::LOTUSginiHybrid<double >;
   getNodePerformance = Helper::purityGini<double >;
   getTermResult = Helper::getMostFreq<double >;
   }


   if (dataSet == 1) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   2,
   //(uli_t)sqrt((double)ncol),
   rng
   );
   } else if (dataSet == 2) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   25,
   //(uli_t)sqrt((double)ncol),
   rng
   );
   } else if (dataSet == 3) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   2,
   //(uli_t)sqrt((double)ncol),
   rng
   );
   } else if (dataSet == 4) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   25,
   //(uli_t)sqrt((double)ncol),
   rng
   );
   } else if (dataSet == 5) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   (uli_t)sqrt((double)ncol),
   rng
   );
   } else if (dataSet == 6) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   2,
   //(uli_t)sqrt((double)ncol),
   rng
   );
   } else if (dataSet == 7) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   25,
   //(uli_t)sqrt((double)ncol),
   rng
   );
   } else if (dataSet == 8) {
   treeCtrl = JTreeCtrl<double >(	getBestVariable,
   getNodePerformance,
   getTermResult,
   (uli_t)sqrt((double)ncol),
   rng
   );
   }

   gettimeofday(&start, 0);
   std::map<uli_t, std::vector<uli_t > > trainDataMap;
   std::vector<uli_t > trainDataRows;
   std::map<uli_t, std::vector<uli_t > > oobDataMap;
   std::vector<uli_t > oobDataRows;

   std::cout 	<< "Grow jungle(" << ntree << "): " << std::endl;
   for(i = 0; i < ntree; ++i) {

   if (dataSet == 2 || dataSet == 4 || dataSet == 8)
   std::cout << i + 1 << std::endl;

   trees.push_back(treeCtrl.makeTree(	data,
   NULL,
   i,
   false,
   //OUTPUT
   trainDataMap,
   trainDataRows,
   oobDataMap,
   oobDataRows
   ));

   if (dataSet == 3)
   std::cout << "tree: " << *trees[i] << std::endl;

   // show freq importance
   treeCtrl.getFreqImportance(trees[i], data, freqMap);

   }
   decGiniMap = treeCtrl.getVarImp();

   gettimeofday(&end, 0);
   std::cout 	<< std::endl << std::endl;
   std::cout 	<< "growing took: " << end.tv_sec-start.tv_sec
   << " sec"
   << std::endl;

   // print freqs
   Helper::inverseMap<uli_t, double>(*freqMap, mapInverse);
   Helper::printTreeMap<double, uli_t, double>(mapInverse, *data, ntree);

   // print importance
   Helper::inverseMap<uli_t, double>(decGiniMap, mapInverse);
   Helper::printTreeMap<double, uli_t, double>(mapInverse, *data, ntree);


   if (dataSet != 3 && dataSet != 7) {
   std::cout << std::endl;
   std::cout << "training set confusion matrix: " << std::endl;
   Helper::showAccur<double >(data, trees, depVar, trainDataMap);

   std::cout << "test set confusion matrix: " << std::endl;
   Helper::showAccur<double >(data, trees, depVar, oobDataMap);
   }

   std::vector<Tree<double, uli_t > * >::iterator itTree(trees.begin());
   while (itTree != trees.end()) {
   delete *itTree;
   ++itTree;
   }

   if (freqMap != NULL) delete freqMap;
   if (data != NULL) delete data;
   */
}

/*
 void TestClass::testRandomJungleCtrl1() const {
 RJungleCtrl<double > rjCtrl;
 gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
 double missingcode = -99;
 uli_t ntree = 1000;
 uli_t dataSet = 10; /////////////
 uli_t depVar = 0;
 uli_t nrow;
 uli_t ncol;
 uli_t varNamesRow = 0;
 uli_t depVarCol = 0;
 uli_t skipRow = 0;
 uli_t skipCol = 0;
 std::string filename;
 std::pair<uli_t, uli_t > csvdim;
 uli_t mtry = 0; // 0 means default value (sqrt)
 GETBESTVARFCT_ASVAR(double);
 GETNODEPERFFCT_ASVAR(double);
 GETTERMRESFCT_ASVAR(double);
 CMPLTREEFCT_ASVAR(double);
 CLASSIFYCMPLDTREEFCT_ASVAR(double);
 getBestVariable = FittingFct::CARTgini<double >;
 //FittingFct::checkUsabilityLOTUSginiHybrid<double >(*data, depVar);
 //getBestVariable = FittingFct::LOTUSginiHybrid<double >;
 getNodePerformance = Helper::purityGini<double >;
 //getTermResult = Helper::getMedianX<double >;
 getTermResult = Helper::getMostFreq<double >;
 cmplTree = CmplFct::T2<double >;
 classifyCmpldTree = CmplFct::classifyT2<double >;

 if (dataSet == 1) {
 filename = "Ger5.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 2;
 } else if (dataSet == 2) {
 filename = "plinkpheno_chrALL_rep0001.csv";
 //filename = "plink_chrALL_rep0001.csv";
 nrow = 3500;//3500;
 ncol = 9188;//9188;
 depVarCol = 5;
 skipRow = 1;
 skipCol = 30;
 mtry = 25;
 //mtry = 0;
 } else if (dataSet == 3) {
 //regression
 filename = "testdatareg.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 2;
 getBestVariable = FittingFct::CARTreg<double >;
 getNodePerformance = Helper::puritySos<double >;
 getTermResult = Helper::getMean2<double >;
 } else if (dataSet == 4) {
 filename = "plink_chrALL_rep0001.csv";
 depVar = 3622;//SNP6_153 0 cM from DR, increases RA risk only in women
 nrow = 3500;//3500;
 ncol = 9100;//9188;
 depVarCol = 5;
 skipRow = 1;
 skipCol = 4;
 mtry = 25;
 }  else if (dataSet == 5) {
 filename = "info.csv";
 //nrow = 1200;
 nrow = 120; // with 120 real samples and } 1000 additional zerosamples
 //, so we avoid the strobel problem
 ncol = 6;    // funny , but true and i dont know why!
 mtry = 2;
 } else if (dataSet == 6) {
 filename = "testdata.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 2;
 } else if (dataSet == 7) {
 filename = "plinkpheno_chrALL_rep0001.csv";
 nrow = 1500;//3500;
 ncol = 9000;//9198;
 depVarCol = 11;
 depVar = 0;//9, = antiCCP (SNP18_269 controls effect of DR on anti-CCP and increases RA risk), 10 = lgM(QTL for SNP11_389),11 = Severty(QTL for SNP9_187/SNP9_192),
 skipRow = 2001;
 skipCol = 20;
 mtry = 0;
 getBestVariable = FittingFct::CARTreg<double >;
 getNodePerformance = Helper::puritySos<double >;
 getTermResult = Helper::getMean2<double >;
 } else if (dataSet == 8) {
 filename = "illumina_50k.csv";
 nrow = 2603;
 ncol = 26267;//26267 only good SNPs //45707 all SNPs;
 depVarCol = 5;
 skipRow = 1;
 skipCol = 6;
 mtry = 0;
 } else if (dataSet == 9) {
 filename = "samples10000.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 10;
 } else if (dataSet == 10) {
 filename = "data.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 }

 DataFrame<double > *data = new DataFrame<double >(nrow, ncol);
 data->setMissingCode(missingcode);
 data->setDepVar(depVar);
 data->getDataFromCSV(filename, nrow, ncol, varNamesRow, depVarCol, skipRow, skipCol);
 if (dataSet == 6) data->printSummary();

 const std::vector<std::string > &varNames = data->getVarNames();
 std::cout << "dep. Var. name:" << varNames[depVar] << std::endl;

 //std::cout << *data << std::endl;

 std::cout 	<< "Grow jungle(" << ntree << "): " << std::endl;

 rjCtrl.makeForest(
 rng,
 getBestVariable,
 getNodePerformance,
 getTermResult,
 cmplTree,
 classifyCmpldTree,
 *data,
 mtry,
 ntree
 );

 std::cout << std::endl << std::endl;

 //rjCtrl.printVarFreqs(*data);
 rjCtrl.printIntrinsicVarImp(*data);
 if (dataSet != 3 && dataSet != 7) rjCtrl.printConfMat(*data);

 if (data != NULL) delete data;
 }

 void TestClass::testRandomJungleCtrl2() const {
 RJungleCtrl<char > rjCtrl;
 gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
 char missingcode = -99;
 uli_t ntree = 5000;
 uli_t dataSet = 8; /////////////
 uli_t depVar = 0;
 uli_t nrow;
 uli_t ncol;
 uli_t varNamesRow = 0;
 uli_t depVarCol = 0;
 uli_t skipRow = 0;
 uli_t skipCol = 0;
 std::string filename;
 std::pair<uli_t, uli_t > csvdim;
 uli_t mtry = 0; // 0 means default value (sqrt)
 GETBESTVARFCT_ASVAR(char);
 GETNODEPERFFCT_ASVAR(char);
 GETTERMRESFCT_ASVAR(char);
 CMPLTREEFCT_ASVAR(char);
 CLASSIFYCMPLDTREEFCT_ASVAR(char);
 getBestVariable = FittingFct::CARTgini<char>;
 //FittingFct::checkUsabilityLOTUSginiHybrid<char >(*data, depVar);
 //getBestVariable = FittingFct::LOTUSginiHybrid<char >;
 getNodePerformance = Helper::purityGini<char>;
 //getTermResult = Helper::getMedianX<char >;
 getTermResult = Helper::getMostFreq<char>;
 cmplTree = CmplFct::T2<char>;
 classifyCmpldTree = CmplFct::classifyT2<char >;

 if (dataSet == 1) {
 filename = "Ger5.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 2;
 } else if (dataSet == 2) {
 filename = "plinkpheno_chrALL_rep0001.csv";
 //filename = "plink_chrALL_rep0001.csv";
 nrow = 3500;//3500;
 ncol = 9188;//9188;
 depVarCol = 5;
 skipRow = 1;
 skipCol = 30;
 mtry = 25;
 //mtry = 0;
 } else if (dataSet == 3) {
 //regression
 filename = "testdatareg.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 2;
 getBestVariable = FittingFct::CARTreg<char >;
 getNodePerformance = Helper::puritySos<char >;
 getTermResult = Helper::getMean2<char >;
 } else if (dataSet == 4) {
 filename = "plink_chrALL_rep0001.csv";
 depVar = 3622;//SNP6_153 0 cM from DR, increases RA risk only in women
 nrow = 3500;//3500;
 ncol = 9100;//9188;
 depVarCol = 5;
 skipRow = 1;
 skipCol = 4;
 mtry = 25;
 }  else if (dataSet == 5) {
 filename = "info.csv";
 //nrow = 1200;
 nrow = 120; // with 120 real samples and } 1000 additional zerosamples
 //, so we avoid the strobel problem
 ncol = 6;    // funny , but true and i dont know why!
 mtry = 2;
 } else if (dataSet == 6) {
 filename = "testdata.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 2;
 } else if (dataSet == 7) {
 filename = "plinkpheno_chrALL_rep0001.csv";
 nrow = 1500;//3500;
 ncol = 9000;//9198;
 depVarCol = 11;
 depVar = 0;//9, = antiCCP (SNP18_269 controls effect of DR on anti-CCP and increases RA risk), 10 = lgM(QTL for SNP11_389),11 = Severty(QTL for SNP9_187/SNP9_192),
 skipRow = 2001;
 skipCol = 20;
 mtry = 0;
 getBestVariable = FittingFct::CARTreg<char >;
 getNodePerformance = Helper::puritySos<char >;
 getTermResult = Helper::getMean2<char >;
 } else if (dataSet == 8) {
 filename = "illumina_50k.csv";
 nrow = 2603;
 ncol = 26267;//26267 only good SNPs //45707 all SNPs;
 depVarCol = 5;
 skipRow = 1;
 skipCol = 6;
 mtry = 0;
 } else if (dataSet == 9) {
 filename = "samples100000.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 mtry = 10;
 } else if (dataSet == 10) {
 filename = "data.csv";
 csvdim = Helper::getCSVdim(filename);
 nrow = csvdim.first;
 ncol = csvdim.second;
 }


 DataFrame<char > *data = new DataFrame<char >(nrow, ncol);
 data->setMissingCode(missingcode);
 data->setDepVar(depVar);
 data->getDataFromCSV(filename, nrow, ncol, varNamesRow, depVarCol, skipRow, skipCol);
 if (dataSet == 6) data->printSummary();

 const std::vector<std::string > &varNames = data->getVarNames();
 std::cout << "dep. Var. name:" << varNames[depVar] << std::endl;

 //std::cout << *data << std::endl;

 std::cout 	<< "Grow jungle(" << ntree << "): " << std::endl;

 rjCtrl.makeForest(
 rng,
 getBestVariable,
 getNodePerformance,
 getTermResult,
 cmplTree,
 classifyCmpldTree,
 *data,
 mtry,
 ntree
 );

 std::cout << std::endl << std::endl;

 //rjCtrl.printVarFreqs(*data);
 rjCtrl.printIntrinsicVarImp(*data);
 if (dataSet != 3 && dataSet != 7) rjCtrl.printConfMat(*data);

 if (data != NULL) delete data;
 }
 */
