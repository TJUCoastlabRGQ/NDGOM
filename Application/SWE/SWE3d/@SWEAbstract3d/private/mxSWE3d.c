#include "mxSWE3d.h"
#include "mex.h"

/** convert mex variable to PhysVolField structure */
PhysField convertMexToPhysField(const mxArray *mxfield) {
    const mwSize *dims = mxGetDimensions(mxfield);
    PhysField field;
    field.Np = dims[0];
    field.K = dims[1];
    field.Nfield = dims[2];
    const size_t Ntmp = field.Np * field.K;
    
    field.h = mxGetPr(mxfield);
    field.hu = field.h + Ntmp;
    field.hv = field.hu + Ntmp;
    field.z = field.hv + Ntmp;
    return field;
}

/** Evaluate the flow rate depending on the depth threshold */
void evaluateFlowRateByDeptheThreshold(
        const double hcrit,  ///< depth threshold
        const double h,      ///< depth
        const double hu,     ///< water flux
        const double hv,     ///< water flux
        double *u,           ///< result velocity
        double *v            ///< velocity
        ) {
    if (h > hcrit) {
        //     const double sqrt2 = 1.414213562373095;
        //     double h4 = pow(h, 4);
        //     *u = sqrt2 * h * hu / sqrt( h4 + max( hcrit, h4 ) );
        //     *v = sqrt2 * h * hv / sqrt( h4 + max( hcrit, h4 ) );
        *u = hu / h;
        *v = hv / h;
    } else {
        *u = 0.0;
        *v = 0.0;
    }
    
    return;
}

/* Evaluate whether the cell adjacent to a given cell is dry*/
void evaluateWetDryInterface(
        signed char *status,
        const mxArray *FToE,
        double *DryFaceFlag
        ){
    const mwSize *dims = mxGetDimensions(FToE);
    double *FaceToElement = mxGetPr(FToE);
    const size_t Ne = dims[1];
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (mwIndex i = 0; i<Ne; i++){
        mwIndex Local_Element = (int)FaceToElement[2 * i];
        mwIndex Adjacent_Element = (int)FaceToElement[2 * i + 1];
        NdgRegionType Local_type = (NdgRegionType)status[Local_Element - 1];
        NdgRegionType Adjacent_type = (NdgRegionType)status[Adjacent_Element - 1];
        if (Local_type != NdgRegionWet || Adjacent_type != NdgRegionWet)
            DryFaceFlag[i] = 1;
    }
}