#ifndef __mxSWE_H__
#define __mxSWE_H__

#include "mex.h"

#if !defined(_WIN32)
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)
#endif

#define EPS 1e-6

typedef enum {
  NdgRegionNormal = 1,
  NdgRegionRefine = 2,
  NdgRegionSponge = 3,
  NdgRegionWet = 4,
  NdgRegionDry = 5,
  NdgRegionPartialWet = 6,
  NdgRegionPartialWetFlood = 7,
  NdgRegionPartialWetDamBreak = 8
} NdgRegionType;

typedef struct {
  size_t Np;      ///< length of 1st dimension
  size_t K;       ///< length of 2nd dimension
  size_t Nfield;  ///< length of 3rd dimension
  double *h;
  double *hu;
  double *hv;
  double *z;
} PhysField;

/** convert mex variable to PhysVolField structure */
 PhysField convertMexToPhysField(const mxArray *);

/** Evaluate the flow rate depending on the depth threshold */
 void evaluateFlowRateByDeptheThreshold(
    const double ,  ///< depth threshold
    const double ,      ///< depth
    const double ,     ///< water flux
    const double ,     ///< water flux
    double *,           ///< result velocity
    double *            ///< velocity
) ;

void evaluateWetDryInterface( signed char *, const mxArray *, double * );

#endif  //__mxSWE_H__
