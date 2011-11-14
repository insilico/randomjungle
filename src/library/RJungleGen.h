#ifndef RJUNGLEGEN_H_
#define RJUNGLEGEN_H_

/*
 * Includes
 */

#ifndef __WINDOWS__
#ifndef __NOPLUGIN__
#include <dlfcn.h> // Guess Standard POSIX/UNIX API
#endif
#else
#ifndef __NOPLUGIN__
#include <windows.h>
#endif
#endif

#include "CmplFct.h"
#include "BuildinGenerators.h"
#include "DataFrame.h"
#include "RJunglePar.h"
#include "ClassAtom.h"
#include "CmpldTree.h"
#include "DataTreeSet.h"
#include "Tree.h"


/*
 * Def.
 */

template <class T >
class RJungleGen{
public:
  RJungleGen() {
  	// call macros to declare important generator functions
    isPlugin = false;
  	wasInitialized = false;
  }

  ~RJungleGen() {
    if (wasInitialized) {
      #ifndef __NOPLUGIN__
      #ifndef STATIC
      #ifndef __WINDOWS__
      if (isPlugin) dlclose(sdl_library);
      #else
      if (isPlugin) FreeLibrary(sdl_library);
      #endif
      #endif
      #endif
    }
  }

  void init(RJunglePar &par, DataFrame<T > &data) {
    wasInitialized = true;

    // define building functions
    if (strcmp(par.plugin, "") != 0) {
      #ifndef __NOPLUGIN__
      isPlugin = true;
      typedef void (*func_t)(ASSIGNMENT_PARAMETERS_PLUGIN(T));
      func_t fptr;

      // check for plugin
      // open plugin
      //verboseStream << "loading plugin..." << std::endl;

      #ifndef STATIC
      #ifndef __WINDOWS__
      sdl_library = dlopen(par.plugin, RTLD_NOW);
      if (sdl_library == NULL) throw(Exception(ERRORCODE_29));

      fptr = (func_t)dlsym(sdl_library, "rjungle_plugin_assignFunctions");
      if (fptr == NULL) throw(Exception(ERRORCODE_30));
      #else
      sdl_library = LoadLibrary(par.plugin);
      if (sdl_library == NULL) throw(Exception(ERRORCODE_29));

      fptr = (func_t)GetProcAddress(sdl_library, "rjungle_plugin_assignFunctions");
      if (fptr == NULL) throw(Exception(ERRORCODE_30));
      #endif
      #endif
      (*fptr)(&this->fct);
      #else
      throw(Exception(ERRORCODE_40));
      #endif
    } else if (par.treeType == tt_CART) {
      rjungle::Generator::
      CART::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CARTreg) {
      rjungle::Generator::
      CARTreg::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CARTregTwoing) {
      rjungle::Generator::
      CARTregTwoing::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_CARTregWithMiss) {
      rjungle::Generator::
      CARTregWithMiss::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CARTtwoing) {
      rjungle::Generator::
      CARTtwoing::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_LOTUS) {
      rjungle::Generator::
      LOTUS::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_UCART) {
      rjungle::Generator::
      UCART::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_CARTsrt) {
      rjungle::Generator::
      CARTsrt::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_SCART) {
      rjungle::Generator::
      SCART::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_trendCART) {
      rjungle::Generator::
      trendCART::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_imputeTree) {
      rjungle::Generator::
      imputeTree::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_twoWayIntTree) {
      rjungle::Generator::
      CARTtwoWayIntAdd::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_CARTcor) {
      rjungle::Generator::
      CARTcor::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CISNP) {
      rjungle::Generator::
      CISNP::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CICase1) {
      rjungle::Generator::
      CICase1::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CICase3) {
      rjungle::Generator::
      CICase3::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CARTsse) { //sse
      rjungle::Generator::
      CARTsse::assignFunctions<T >(&this->fct);
    } else {
      rjungle::Generator::
      CART::assignFunctions<T >(&this->fct);
  	}

    // check usability of trees regarding data
    if (fct.checkUsability != NULL) (*fct.checkUsability)(data);
  }

  BuildinGenFct<T > fct;

  // for plugins
  #ifndef __NOPLUGIN__
  #ifndef __WINDOWS__
  void* sdl_library;
  #else
  HMODULE sdl_library;
  #endif
  #endif

  bool isPlugin;
  bool wasInitialized;
};

#endif /* RJUNGLEGEN_H_ */
