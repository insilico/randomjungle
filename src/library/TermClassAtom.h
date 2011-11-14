#ifndef TERMCLASSATOM_H_
#define TERMCLASSATOM_H_


#include <vector>
#include <limits>
#include "ClassAtom.h"
#include "treedefs.h"

template <class T, class C >
class TermClassAtom : public ClassAtom<T, C > {
public:
	TermClassAtom(): val(T()) { };
	TermClassAtom(const TermClassAtom<T, C > &termCA) : val(termCA.val) { };
	TermClassAtom(T val) : val(val) { };

	virtual ~TermClassAtom(void) {
	};

	virtual inline T getResult() const {
		return this->val;
	}

	virtual inline bool isConsistent() const {
		return true;
	}



	/*
	 */
	virtual std::ostream& print(std::ostream& os) const {
		os << "term:" << (double)this->val;
		return os;
	}

	virtual void printXml(std::ostream& os) const {
    os
    << "<classifier id=\"0\" size=\"1\" type=\"term\">"
    << "<value id=\"0\" type=\"eq\">" << (double)this->val << "</value>"
		<< "</classifier>" << std::endl;
	}

	virtual inline uli_t getType() const { return ca_TERM; }

	T val;

};

#endif /*TERMCLASSATOM_H_*/
