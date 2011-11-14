#ifndef LOTUSTERMCLASSATOM_H_
#define LOTUSTERMCLASSATOM_H_


#include <vector>
#include <limits>
#include "TermClassAtom.h"
#include "treedefs.h"

template <class T, class C >
class LotusTermClassAtom : public TermClassAtom<T, C > {
public:
	LotusTermClassAtom(): val(T()) { };
	LotusTermClassAtom(const LotusTermClassAtom<T, C > &termCA) : val(termCA.val) { };
	LotusTermClassAtom(T val) : val(val) { };

	virtual ~LotusTermClassAtom(void) {
	};

	/// waste
	T val;

	/// x predictor variable
	uli_t varID;

	/// Betas of logistic regression
  std::vector<double > betas;
};

#endif /*LOTUSTERMCLASSATOM_H_*/
