#ifndef MXSWE1D_H
#define MXSWE1D_H

#include "mex.h"
#include <math.h>

#define TAU 1e-6

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define EPS 1e-6
#define NVAR 2

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

typedef enum {
  NdgEdgeInner = 0,
  NdgEdgeGaussEdge = 1,
  NdgEdgeSlipWall = 2,
  NdgEdgeNonSlipWall = 3,
  NdgEdgeZeroGrad = 4,
  NdgEdgeClamped = 5,
  NdgEdgeClampedDepth = 6,
  NdgEdgeClampedVel = 7,
  NdgEdgeFlather = 8
} NdgEdgeType;

typedef struct {
  size_t Np;     ///< length of 1st dimension
  size_t K;      ///< length of 2nd dimension
  size_t Nfield; ///< length of 3rd dimension
  double *h;
  double *hu;
  double *z;
} PhysField1d;

/** convert mex variable to PhysVolField structure */
PhysField1d convertMexToPhysField(const mxArray * ); 

/** Evaluate the flow rate depending on the depth threshold */
void evaluateFlowRateByDeptheThreshold(const double , ///< depth threshold
                                  const double ,     ///< depth
                                  const double ,    ///< water flux
                                  double *           ///< result velocity
);

void evaluateWetDryInterface( signed char *, const mxArray *, double * );

#endif // MXSWE1D_H