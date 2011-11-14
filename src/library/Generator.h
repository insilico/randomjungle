#ifndef GENERATOR_H_
#define GENERATOR_H_

/*
 * Includes
 */

#include <iostream>
#include <cstring>
#include <vector>
#include <string>
#include <limits>
#include <ctime>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Helper.h"
#include "CmplFct.h"
#include "treedefs.h"


/*
 * Def.
 */

template <class T >
class Generator {
public:
  // const. / dest.
	Generator() { }
	virtual ~Generator() { }
	
	// find the best variable in each node
	virtual static GETBESTVARFCT(T) {

  }
  
	virtual static GETNODEPERFFCT(T) {
	
  }
  
	virtual static GETTERMRESFCT(T) {
	
	}
	
	virtual static CMPLTREEFCT(T) {
	
	}

  virtual static CLASSIFYCMPLDTREEFCT(T) {
  
  }
}

// definition for plugins
// not supported yet
/*
#define PLUGIN_INIT(class_name)\
	extern "C"\
	Generator* _librjungle_plugin_new() {\
		return new class_name();\
	}\
	extern "C"\
	void _librjungle_plugin_delete(Generator* plugin) {\
		delete plugin;\
	}
*/
#endif /* GENERATOR_H */
