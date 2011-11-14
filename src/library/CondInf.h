/* 
 * Copyright (C) 2008-2010  Daniel F. Schwarz
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

#ifndef CONDINF_H_
#define CONDINF_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include "Helper.h"
#include "FittingFct.h"
#include "mvt.h"

/** 
 * \brief CI Functions for various cases as follows:
 * - y (numerical), x(numerical)
 * - y (nominal), x(numerical)
 * - SNP data: y (2-classes), x(numerical)
 */
template<class T>
class CondInf {
public:
  CondInf():
    xIsNominal(false), yIsNominal(false), p(1), q(1), maxpts(25000),
      tol(1e-10), abseps(0.001), releps(0.0), error(0.0), rng(NULL) {}

  virtual ~CondInf(){};

  /*
   * Conditional Inference Trees for classification
   * ----------------------------------------------
   *
   * Variables:
   *   y   \in (nominal)
   *   x_j \in (numeric/ordinal)
   *
   */
  ClassAtom<T, uli_t> *CIClassification(FittingFctPar<T> &fitFctPar) {

    // fitting params
    std::vector<uli_t> &colMask = *fitFctPar.colMaskVec;
    std::vector<uli_t> &rowMask = *fitFctPar.rowMaskVec;
    std::vector<double> &rowWeight = *fitFctPar.rowWeight;
    DataFrame<T> &data = *fitFctPar.data;
    uli_t depVar = fitFctPar.depVar;

    // variables
    T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
    double bestPValue = 2.0;
    double alpha = 0.05;
    uli_t i, bestVar, curVar;
    bestVar = curVar = 0;
    rng = fitFctPar.rng;

    /*
     * Assign plugin parameters
     */
    std::vector<std::string> pluginPar;
    std::string strDump(data.par.pluginPar);
    Helper::strSplit(strDump, ",", pluginPar);

    if (pluginPar.size() > 0)
      if (pluginPar[0] != "")
        alpha = atof(pluginPar[0].c_str());

    /*
     * Optional feature: First node splits always variable xyz
     */
    bool isCorrectedTree = false;
    bool isFirstVar = (fitFctPar.depth == 0) ? true : false;
    uli_t idxFirstVar = 0;

    // Load 1st var.
    if (isFirstVar)
      if (pluginPar.size() > 1)
        if (pluginPar[1] != "") {
          idxFirstVar = data.getIndexOfVarName(pluginPar[1]);
          isCorrectedTree = true;
        }

    /*
     * init 
     */
    data.getRowMaskedCol(depVar, rowMask, y);
    Helper::getMaskedVec(rowWeight, rowMask, w);

    prepareY(*this, data, data.par.depVar);

    /* if only one class left, stop */
    if (yIsNominal && (q < 2)) {
      fitFctPar.stopGrowing = true;
      return classifier;
    }

    /*
     * global test
     */
    if (isFirstVar && isCorrectedTree) {
      curVar = idxFirstVar;
      bestVar = idxFirstVar;
    } else {
      for (i = 0; i < colMask.size(); ++i) {
        //skip dep. var.
        if (colMask[i] == depVar) {
          continue;
        }

        curVar = colMask[i];

        // get x
        data.getRowMaskedCol(curVar, rowMask, x);

        prepareX(*this, data, curVar);

        // calc. p value
        getPAsym(*this);

        if (xIsNominal && (p < 2))
          continue;

        if (pAsym < bestPValue) {
          bestPValue = pAsym;
          bestVar = curVar;
        }
      }

      // check global hypo. (bonferroni corrected)
      //if (bestPValue > alpha / colMask.size()) {
      if (bestPValue > alpha) {
        fitFctPar.stopGrowing = true;
        return classifier;
      }
    }

    // get split
    data.getRowMaskedCol(bestVar, rowMask, x);
    getSplit(*this);

    // build classifier
    classifier->setVarID(bestVar);
    classifier->setThreshold(split);
    classifier->setMissingCode(data.getMissingCode());

    *fitFctPar.treeImprovement = 1;
    *fitFctPar.bestVar = bestVar;
    return classifier;
  }

  /*
   * Conditional Inference Trees for
   * -------------------------------
   *
   * Variables:
   *   y   \in (nominal) -> h(y,...) = e_k(y)
   *   x_j \in (numeric/ordered) -> g(x) = x or rank(x)
   *
   * Global hypo.:
   *   - Using only "standardized linear statistic":
   *       c_max = max|(t-\mu)_k / \sqrt(\dCov_T_{kk})|
   *   - Bonferroni correction for mult. testing (unused)
   *
   * Splitting:
   *   - Find best split point in ordered subset A. So, find split in O(n).
   */
  static ClassAtom<T, uli_t> *CICase1(FittingFctPar<T> &fitFctPar) {
		CondInf<T> ci;
    ci.prepareY = prepareYCase1;
    ci.prepareX = prepareXCase1;
    ci.getPAsym = getPAsymCase1;
    ci.getSplit = calcSplitCase1;
    return ci.CIClassification(fitFctPar);
  }

  static void prepareYCase1(CondInf<T> &ci, DataFrame<T> &data, uli_t j) {
    ci.calcSumw();
    ci.trafoRank(ci.y);
    ci.yIsNominal = true;
    ci.setQ();
    ci.calcdExp_y();
    ci.calcdCov_y();
  }

  static void prepareXCase1(CondInf<T> &ci, DataFrame<T> &data, uli_t j) {
    ci.xIsNominal = false;
    ci.setP();
    if (data.lom[j] == lom_ordinal)
      ci.trafoRank(ci.x);
    ci.calcSumwxCT1();
    ci.calcdCov_T();
    ci.calcdExp_T();
    ci.calcT();
  }

  static void getPAsymCase1(CondInf<T> &ci) {
    ci.calcCMax();
    ci.C_maxabsConditionalPvalue();
  }

  static void checkForCase1(DataFrame<T> &data) {
    checkForCI(data);
  }

  static void calcSplitCase1(CondInf<T> &ci) {
    std::vector<T> levels, xBackup;
    double cMax = -1;
    uli_t i, j;

    // initial calc.
    ci.calcSumw();
    ci.trafoRank(ci.y);
    ci.yIsNominal = true;
    ci.setQ();
    ci.calcdExp_y();
    ci.calcdCov_y();

    // get unique values
    Helper::getClasses(ci.x, levels);
    sort(levels.begin(), levels.end());

    // max. lin. statistic
    xBackup = ci.x;
    for (i = 0; i < levels.size() - 1; ++i) {

      // prepare x
      for (j = 0; j < ci.x.size(); ++j) {
        ci.x[j] = (xBackup[j] <= levels[i]) ? 1 : 0;
      }

      // calc. cMax
      ci.calcSumwxCT1();
      ci.calcdCov_T();
      ci.calcdExp_T();
      ci.calcT();
      ci.calcCMax();

      if (ci.cMax > cMax) {
        cMax = ci.cMax;
        ci.split = levels[i];
      }
    }
    ci.x = xBackup;
  }

  /*
   * calc. asym. p value for case 1: y (nominal) / x (numeric/ordinal)
   */
  void calcPAsymCase1() {
    calcSumw();
    calcSumwxCT1();
    trafoRank(y);
    yIsNominal = true;
    setQ();
    calcdExp_y();
    calcdCov_y();
    calcdCov_T();
    calcdExp_T();
    calcT();

    calcCMax();
    calcPAsym();
  }

  /**
	 * \brief
   * Conditional Inference Trees for SNPs
   * ------------------------------------
   *
   * Variables:
   *   y   \in {0,1}   (nominal) -> h(y,...) = e_2(y)
   *   x_j \in {0,1,2} (numeric) ->     g(x) = x
   *
   * Global hypo.:
   *   - Using only "standardized linear statistic":
   *       c_max = max|(t-\mu)_k / \sqrt(\dCov_T_{kk})|
   *   - Bonferroni correction for mult. testing (unused)
   *
   * Splitting:
   *   - Find best split point in ordered subset A. So, find split in O(n).
   */
  static ClassAtom<T, uli_t> *CISNP(FittingFctPar<T> &fitFctPar) {
		CondInf<T> ci;
    ci.prepareY = prepareYCase2;
    ci.prepareX = prepareXCase2;
    ci.getPAsym = getPAsymCase2;
    ci.getSplit = calcSplitCase2;
    return ci.CIClassification(fitFctPar);
  }

	static void prepareYCase2(CondInf<T> &ci, DataFrame<T> &data, uli_t j) {
    ci.calcSumw();
    ci.trafoRank(ci.y);
    ci.yIsNominal = true;
    ci.setQ();
  }

	static void prepareXCase2(CondInf<T> &ci, DataFrame<T> &data, uli_t j) {
    ci.xIsNominal = false;
    ci.setP();
    if (data.lom[j] == lom_ordinal)
      ci.trafoRank(ci.x);
  }

	static void getPAsymCase2(CondInf<T> &ci) {
    ci.calcCMaxCase2();
    ci.calcPAsymCase2();
  }

  static void checkForTwoSampleData(DataFrame<T> &data) {
    if (data.varCategories[data.par.depVar].size() != 2)
      throw Exception(ERRORCODE_48);

    checkForCI(data);
  }

	static void calcSplitCase2(CondInf<T> &ci) {
    std::vector<T> levels, xBackup;
    double cMax = -1;
    uli_t i, j;

    // initial calc.
    ci.calcSumw();
    ci.trafoRank(ci.y);
    ci.yIsNominal = true;
    ci.setQ();

    // get unique values
    Helper::getClasses(ci.x, levels);
    sort(levels.begin(), levels.end());

    // max. lin. statistic
    xBackup = ci.x;
    for (i = 0; i < levels.size() - 1; ++i) {

      // prepare x
      for (j = 0; j < ci.x.size(); ++j) {
        ci.x[j] = (xBackup[j] <= levels[i]) ? 1 : 0;
      }

      // calc. cMax
      ci.calcCMaxCase2();

      if (ci.cMax > cMax) {
        cMax = ci.cMax;
        ci.split = levels[i];
      }
    }
    ci.x = xBackup;
  }

  void calcCMaxCase2() {
    double wx1, wx2, wx, w1, w2, sumw, w_star, s_star, sumwxsqr, cmax;
    wx1 = wx2 = wx = w1 = w2 = sumw = w_star = s_star = sumwxsqr = cmax = 0;

    for (uli_t i = 0; i < x.size(); ++i) {
      if (y[i] == y[0]) {
        wx1 += w[i] * x[i];
        w1 += w[i];
      } else {
        wx2 += w[i] * x[i];
        w2 += w[i];
      }

      sumwxsqr += pow(w[i] * x[i], 2);
    }

    wx = wx1 + wx2;
    sumw = w1 + w2;

    w_star = w1 / sumw;

    s_star = (w_star - pow(w_star, 2)) / (sumw - 1);
    s_star *= (sumw * sumwxsqr - pow(wx, 2));

    cMax = fabs((wx1 - w_star * wx) / sqrt(s_star));
  }

  void calcPAsymCase2() {
    pAsym = gsl_cdf_ugaussian_P(-1.0 * fabs(cMax)) * 2.0;
  }


  /**
	 * \brief
   * Conditional Inference Trees for
   * -------------------------------
   *
   * Variables:
   *   y   \in (numeric) -> h(y,...) = y
   *   x_j \in (numeric) -> g(x) = x
   *
   * Global hypo.:
   *   - Using only "standardized linear statistic":
   *       c_max = max|(t-\mu)_k / \sqrt(\dCov_T_{kk})|
   *   - Bonferroni correction for mult. testing (unused)
   *
   * Splitting:
   *   - Find best split point in ordered subset A. So, find split in O(n).
   */
  static ClassAtom<T, uli_t> *CICase3(FittingFctPar<T> &fitFctPar) {
		CondInf<T> ci;
    ci.prepareY = prepareYCase3;
    ci.prepareX = prepareXCase3;
    ci.getPAsym = getPAsymCase3;
    ci.getSplit = calcSplitCase3;
    return ci.CIClassification(fitFctPar);
  }

	static void prepareYCase3(CondInf<T> &ci, DataFrame<T> &data, uli_t j) {
    ci.calcSumw();
    ci.yIsNominal = false;
    ci.setQ();
    ci.calcdExp_y();
    ci.calcdCov_y();
  }

	static void prepareXCase3(CondInf<T> &ci, DataFrame<T> &data, uli_t j) {
    ci.xIsNominal = false;
    ci.setP();
    if (data.lom[j] == lom_ordinal)
      ci.trafoRank(ci.x);
    ci.calcSumwxCT1();
    ci.calcdCov_T();
    ci.calcdExp_T();
    ci.calcT();
  }

	static void getPAsymCase3(CondInf<T> &ci) {
    ci.calcCMax();
    ci.C_maxabsConditionalPvalue();
  }

  static void checkForCase3(DataFrame<T> &data) {
    checkForCI(data);
  }

	static void calcSplitCase3(CondInf<T> &ci) {
    std::vector<T> levels, xBackup;
    double cMax = -1;
    uli_t i, j;

    // initial calc.
    ci.calcSumw();
    ci.yIsNominal = false;
    ci.setQ();
    ci.calcdExp_y();
    ci.calcdCov_y();

    // get unique values
    Helper::getClasses(ci.x, levels);
    sort(levels.begin(), levels.end());

    // max. lin. statistic
    xBackup = ci.x;
    for (i = 0; i < levels.size() - 1; ++i) {

      // prepare x
      for (j = 0; j < ci.x.size(); ++j) {
        ci.x[j] = (xBackup[j] <= levels[i]) ? 1 : 0;
      }

      // calc. cMax
      ci.calcSumwxCT1();
      ci.calcdCov_T();
      ci.calcdExp_T();
      ci.calcT();
      ci.calcCMax();

      if (ci.cMax > cMax) {
        cMax = ci.cMax;
        ci.split = levels[i];
      }
    }
    ci.x = xBackup;
  }

  /*
   *
   * Helper functions
   * for all cases
   *
   */

  // get sum of w
  void calcSumw() {
    sumw = 0;
    for (uli_t i = 0; i < w.size(); ++i) {
      sumw += w[i];
    }
  }

  // get sum of w * x
  void calcSumwxCT1() {

    assert(w.size() == x.size());

    swx.clear();
    swx.assign(p, 0);

    CT1.clear();
    CT1.assign(p * p, 0);

    uli_t i;
    double tmp;

    if (xIsNominal) {
      // x is nominal
      for (i = 0; i < x.size(); ++i) {
        if (w[i] == 0)
          continue;

        swx[(size_t)x[i]] += w[i];
        CT1[(size_t)(x[i] * p + x[i])] += w[i];
      }
    } else {
      // x is numerical
      assert(p == 1);

      for (i = 0; i < x.size(); ++i) {
        if (w[i] == 0)
          continue;

        tmp = w[i] * x[i];
        swx[0] += tmp;
        CT1[0] += tmp * x[i];
      }
    }

  }

  // trafo. values to rank
  void trafoRank(std::vector<T> &vec) {
    std::vector<std::pair<T, uli_t> > freqPairs;
    uli_t i;

    for (i = 0; i < vec.size(); ++i)
      freqPairs.push_back(std::pair<T, uli_t>(vec[i], i));

    sort(freqPairs.begin(), freqPairs.end());

    if (vec.size() < 1)
      throw Exception("vec.size() < 1");

    T lastVal = freqPairs[0].first;
    T lastRank = 0;
    for (i = 0; i < freqPairs.size(); ++i) {
      if (vec[freqPairs[i].second] != lastVal) {
        lastRank++;
        lastVal = vec[freqPairs[i].second];
      }
      vec[freqPairs[i].second] = lastRank;
    }

  }

  // calc. and set q
  //
  void setQ() {
    uli_t i;

    if (!yIsNominal) {
      q = 1;
    } else {
      T max = 0;
      for (i = 0; i < y.size(); ++i)
        if (max < y[i])
          max = y[i];

      q = (uli_t)(max + 1);
    }
  }

  // calc. and set p
  //
  void setP() {
    uli_t i;

    if (!xIsNominal) {
      p = 1;
    } else {
      T max = 0;
      for (i = 0; i < x.size(); ++i)
        if (max < x[i])
          max = x[i];

      p = (uli_t)(max + 1);
    }
  }

  // get expectation of y
  void calcdExp_y() {
    uli_t i;

    dExp_y.clear();
    dExp_y.assign(q, 0);

    if (yIsNominal) {
      // y is nominal
      for (i = 0; i < y.size(); ++i) {
        if (w[i] == 0)
          continue;
        dExp_y[(size_t)y[i]] += w[i];
      }
    } else {
      // y is numerical
      assert(q == 1);

      for (i = 0; i < y.size(); ++i) {
        if (w[i] == 0)
          continue;
        dExp_y[0] += w[i] * y[i];
      }
    }

    for (i = 0; i < dExp_y.size(); ++i) {
      dExp_y[i] /= sumw;
    }
  }

  // calc variance of y
  void calcdCov_y() {
    uli_t i, j, k, jq;
    double tmp, tmpUnit;

    //calcdExp_y();

    // init. matrix
    dCov_y.clear();
    dCov_y.assign(q * q, 0);

    // calc. variance
    if (yIsNominal) {
      // y is nominal
      for (i = 0; i < y.size(); ++i) {
        if (w[i] == 0)
          continue;

        for (j = 0; j < q; j++) {
          tmpUnit = (y[i] == (T) j) ? 1 : 0;
          tmp = w[i] * (tmpUnit - dExp_y[j]);
          jq = j * q;
          for (k = 0; k < q; k++) {
            tmpUnit = (y[i] == (T) k) ? 1 : 0;
            dCov_y[jq + k] += tmp * (tmpUnit - dExp_y[k]);
          }
        }
      }
    } else {
      // y is numerical
      assert(q == 1);

      for (i = 0; i < y.size(); ++i) {
        if (w[i] == 0)
          continue;

        dCov_y[0] += w[i] * pow(y[i] - dExp_y[0], 2);
      }
    }

    for (i = 0; i < q * q; ++i)
      dCov_y[i] /= sumw;
  }

  // calc variance of y
  void calcdCov_T() {
    double f1, f2;
    uli_t i, pq;
    assert(x.size() == w.size());

    pq = p * q;
    dCov_T.clear();
    dCov_T.assign(pq * pq, 0);
    f2 = 1 / (sumw - 1);
    f1 = sumw * f2;

    if (pq == 1) {
      dCov_T[0] = f1 * dCov_y[0] * CT1[0];
      dCov_T[0] -= f2 * dCov_y[0] * swx[0] * swx[0];
    } else if (p == 1) {
      double swx2, sum2wx, sumw_1;
      swx2 = sum2wx = 0;
      for (i = 0; i < x.size(); ++i) {
        swx2 += pow(x[i] * w[i], 2);
      }
      sum2wx = pow(swx[0], 2);

      // get dCov_y
      assert(dCov_y.size() == q * q);

      // get sum of weights
      sumw_1 = sumw - 1;

      // get dCov_T
      for (i = 0; i < q * q; ++i) {
        dCov_T[i] = dCov_y[i] / sumw_1;
        dCov_T[i] *= sumw * swx2 - sum2wx;
      }
    } else {
      /* two more helpers needed */
      std::vector<double> CT2(pq * pq, 0);
      std::vector<double> Covy_x_swx(pq * q, 0);

      C_kronecker(dCov_y, q, q, CT1, p, p, dCov_T);
      C_kronecker(dCov_y, q, q, swx, p, 1, Covy_x_swx);
      C_kronecker(Covy_x_swx, pq, q, swx, 1, p, CT2);

      for (i = 0; i < (pq * pq); i++)
        dCov_T[i] = f1 * dCov_T[i] - f2 * CT2[i];

    }

  }

  void calcdExp_T() {
    uli_t i, j;

    dExp_T.clear();
    dExp_T.assign(p * q, 0);

    for (i = 0; i < p; ++i)
      for (j = 0; j < q; ++j)
        dExp_T[j * p + i] = swx[i] * dExp_y[j];
  }

  void calcT() {
    t.clear();
    t.assign(p * q, 0);

    if (xIsNominal && yIsNominal)
      for (uli_t i = 0; i < y.size(); ++i)
        t[(size_t)(y[i] * p + x[i])] += w[i];

    if (!xIsNominal && yIsNominal)
      for (uli_t i = 0; i < y.size(); ++i)
        t[(size_t)y[i]] += w[i] * x[i];

    if (xIsNominal && !yIsNominal)
      for (uli_t i = 0; i < y.size(); ++i)
        t[(size_t)x[i]] += w[i] * y[i];

    if (!xIsNominal && !yIsNominal)
      for (uli_t i = 0; i < y.size(); ++i)
        t[0] += w[i] * x[i] * y[i];
  }

  void calcCMax() {
    double c, sd;
    uli_t i, pq;

    c = sd = 0;
    pq = p * q;

    for (i = 0; i < pq; ++i) {
      sd = dCov_T[i * pq + i];
      c = fabs(t[i] - dExp_T[i]) / sqrt(sd);

      if (c > cMax)
        cMax = c;
    }
  }

  void calcPAsym() {
    pAsym = gsl_cdf_ugaussian_P(-cMax) * 2;
  }

  /*
   * Borrowed from Hothorn's party package
   */
  double C_maxabsConditionalPvalue() {

    double tstat = cMax;
    int pq = p * q;
    gsl_rng *r = rng;

    assert(r != NULL);

    int i, j, *infin, sub, n, nu, inform;
    double *lower, *upper, *delta, *corr, *sd, ans, error, prob;

    /* univariate problem */
    if (pq == 1)
      return (gsl_cdf_ugaussian_P(-1.0 * fabs(tstat)) * 2.0);

    /* define variables */
    nu = 0;
    n = pq;

    if (n == 2)
      corr = new double[1];
    else
      corr = new double[n + ((n - 2) * (n - 1)) / 2];

    sd = new double[n];
    lower = new double[n];
    upper = new double[n];
    infin = new int[n];
    delta = new double[n];

    /* mvtdst assumes the unique elements of the triangular
     covariance matrix to be passes as argument CORREL
     */

    for (i = 0; i < n; i++) {
      /* standard deviations */
      if (dCov_T[i * n + i] < tol)
        sd[i] = 0.0;
      else
        sd[i] = sqrt(dCov_T[i * n + i]);

      /* always look at the two-sided problem */
      lower[i] = fabs(tstat) * -1.0;
      upper[i] = fabs(tstat);
      infin[i] = 2;
      delta[i] = 0.0;
      for (j = 0; j < i; j++) {
        sub = (int)floor(j + 1 + (double) ((i - 1) * (i)) / 2 - 1);
        if (sd[i] == 0.0 || sd[j] == 0.0) {
          corr[sub] = 0.0;
        } else {
          corr[sub] = dCov_T[i * n + j] / (sd[i] * sd[j]);
        }
      }
    }

    prob = 0.0;

    /* call FORTRAN subroutine */
    /* Artifacts arises if multiple threads simultaneously call fortran code
     * with sub functions.
     * http://groups.google.com/group/comp.lang.fortran/browse_thread/thread/d755064b5255c587
     */
#pragma omp critical
    {
      mvtdst_(
        &n, &nu, lower, upper, infin, corr, delta, &maxpts, &abseps, &releps,
        &error, &prob, &inform, &r);
    }

    /* inform == 0 means: everything is OK */
    inform = inform;

    switch (inform) {
      case 0:
        break;
      case 1:
        std::cout << "cmvnorm: completion with ERROR > EPS\n";
        break;
      case 2:
        std::cout << "cmvnorm: N > 1000 or N < 1\n";
        prob = 0.0;
        break;
      case 3:
        std::cout << "cmvnorm: correlation matrix not positive semi-definite\n";
        prob = 0.0;
        break;
      default:
        std::cout << "cmvnorm: unknown problem in MVTDST\n";
        prob = 0.0;
    }
    ans = prob;

    /* free your mind ;) */
    delete[] corr;
    delete[] sd;
    delete[] lower;
    delete[] upper;
    delete[] infin;
    delete[] delta;

    pAsym = 1 - ans;

    return (1 - ans); /* return P-value */
  }

  /*
   * Borrowed from Hothorn's party package
   */
  void C_kronecker(
    std::vector<double> &A, uli_t m, uli_t n, std::vector<double> &B, uli_t r,
    uli_t s, std::vector<double> &ans) {

    uli_t i, j, k, l, mr, js, ir;
    double y;

    mr = m * r;
    for (i = 0; i < m; i++) {
      ir = i * r;
      for (j = 0; j < n; j++) {
        js = j * s;
        y = A[j * m + i];
        for (k = 0; k < r; k++) {
          for (l = 0; l < s; l++) {
            ans[(js + l) * mr + ir + k] = y * B[l * r + k];
          }
        }
      }
    }
  }

  static void checkForCI(DataFrame<T> &data) {
    if (data.par.targetPartitionSize < 20)
      throw Exception(ERRORCODE_49);
  }

	/// Modify y vector (generic function)
  void (*prepareY)(CondInf<T> &ci, DataFrame<T> &, uli_t row);
	
	/// Modify y vector (generic function)
  void (*prepareX)(CondInf<T> &ci, DataFrame<T> &, uli_t row);

	/// Calc. asymtotic p-value (generic function)
  void (*getPAsym)(CondInf<T> &ci);
	
	/// Get cutoff/split (generic function)
  void (*getSplit)(CondInf<T> &ci);


	/// Data set is partioned at this value
	T split;

	/// Indicators / level of measurement of variables
  bool xIsNominal, yIsNominal;

	/// Number of classes in x
  uli_t p;

	/// Number of classes in y
	uli_t q;
	
	/// Sum of weights
  double sumw;

	/*!
	 * Value of test statistics c.
	 * Here, it is the maximum of the absolute values of the standardized linear statistic.j
	 */
	double cMax;

	/// P-value of test 
	double pAsym;

	/// Vectors of x values
  std::vector<T> x;

	/// Vector of y values
	std::vector<T> y;

  /// Weights
  std::vector<double> w;

	/// Measurement of association between y and x via multivariate linear statistic t.
	std::vector<double> t;
	
	/// Conditional expectation of T statistic.
	std::vector<double> dExp_T;

  /// Covariance matrix of T statistic.
	std::vector<double> dCov_T;

	/// Conditional expectation of the influence function
	std::vector<double> dExp_y;
	
	/// Conditional covariance of the influence function
	std::vector<double> dCov_y;

	/// Sum of the case weights * x
	std::vector<double> swx;

  /// Helping structure
	std::vector<double> CT1;

	/// Parameter of Fortran function mvtdst
  int maxpts;

	/// Parameter of Fortran function mvtdst
	int inform;

	/// Parameter of Fortran function mvtdst
  double tol;

	/// Parameter of Fortran function mvtdst
	double abseps;

	/// Parameter of Fortran function mvtdst
	double releps;

	/// Parameter of Fortran function mvtdst
	double error;

	/// Random number generator
  gsl_rng *rng;
};

#endif /* CONDINF_H_ */
