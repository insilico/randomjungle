/*
 * JTreeTweaks.h
 *
 *  Created on: 08.02.2009
 *      Author: schwarz
 */

#ifndef JTREETWEAKS_H_
#define JTREETWEAKS_H_

typedef unsigned long int uli_t;

#ifdef __cplusplus
extern "C" {
#endif

  // UNUSED STRUCT ! SEE RJunglePar.h
  typedef struct JTreeTweaks {
    uli_t maxTreeDepth;
    uli_t targetPartitionSize;
  } JTreeTweaks;

#ifdef __cplusplus
}
#endif

#endif /* JTREETWEAKS_H_ */
