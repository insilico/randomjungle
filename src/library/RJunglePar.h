#ifndef RANDOMJUNGLEPAR_H_
#define RANDOMJUNGLEPAR_H_

#include <gsl/gsl_rng.h>
#include "treedefs.h"

/*
 * Def.
 */
#ifndef HAVE__BOOL
#define bool int
#endif

//typedef unsigned long int uli_t;

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct RJunglePar {
    char *filename;
    char delimiter;
    unsigned int treeType;
    uli_t ntree;
    uli_t ntreeMpi; // Total number of trees when MPI is used
    uli_t mtry;
    unsigned int depVar;
    char *depVarName;
    uli_t nrow;
    uli_t ncol;
    unsigned int varNamesRow;
    unsigned int depVarCol;
    char *outprefix;
    uli_t skipRow;
    uli_t skipCol;
    int missingcode;
    unsigned int impMeasure;
    unsigned int backSel;
    uli_t numOfImpVar;
    bool downsampling_flag;
    bool verbose_flag;
    unsigned int memMode;
    unsigned int saveJungleType;
    char *predict;
    int varproximities;
    bool summary_flag;
    bool testlib_flag;
    char *plugin;
    char *colSelection;
    unsigned int imputeIt;
    bool gwa_flag;
    bool allcont_flag;
    bool transpose_flag;
    bool sampleproximities_flag;
    bool weightsim_flag;
    bool extractdata_flag;
    bool yaimp_flag;
    char *classweights;
    unsigned int seed;
    int nthreads;
    bool pedfile_flag;
    char *pluginPar;
    double tunemtry;
    int outlier;
    bool votes_flag;
    bool oob_flag;
    int prototypes;
    int mpi;
    int mpiId;
    int mpiSize;
    double condimp;
    bool permresponse_flag;

    char delimScale;
    double cutoffHighLD;

    gsl_rng *rng;
    char *version;

    // JTree tweaks
    uli_t maxTreeDepth;
    uli_t targetPartitionSize;
  } RJunglePar;

#ifdef __cplusplus
}
#endif

#endif /*RANDOMJUNGLEPAR_H_*/
