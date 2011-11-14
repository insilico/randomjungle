#ifndef LIBRANDOMJUNGLE_H
#define LIBRANDOMJUNGLE_H

/*
 * Includes
 */


#include "treedefs.h"
#include "RJunglePar.h"
#include "TimeProf.h"

/*
 * Preprocessor commands / marcos
 */

#ifdef __cplusplus
extern "C" {
#endif

//__lrj void __cdecl randomJungle(RJunglePar);
//__lrj RJunglePar __cdecl initRJunglePar();
//__lrj void __cdecl testLibrandomjungle();

void __cdecl randomJungle(RJunglePar);
RJunglePar __cdecl initRJunglePar();
void __cdecl testLibrandomjungle();


#ifdef __cplusplus
}
#endif

#endif  // LIBRANDOMJUNGLE_H
