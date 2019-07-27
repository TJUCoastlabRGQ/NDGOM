#ifndef __mxNonhydro_SWE_H__
#define __mxNonhydro_SWE_H__

typedef enum {
    Inner = 0,
    GaussEdge = 1,
    SlipWall = 2,
    NonSlipWall = 3,
    ZeroGrad = 4,
    Clamped = 5,
    ClampedDepth = 6,
    ClampedVel = 7,
    Flather = 8,
    NonLinearFlatherDepth = 9,
    NonLinearFlatherFlow = 10,
    NonReflectFlux = 11,
    BottomBoundary  = 12,
    UpperSurfaceBoundary =13
} NdgEdgeType;

void mxSetEdgeType( NdgEdgeType EdgeType, double * EdgeIndex ){
        switch(EdgeType){
            /*Face type at the wall boundary is set to be the Neumann boundary, and is flagged by 1*/
            case GaussEdge:
                
            case SlipWall:
                *EdgeIndex = 1;
                break;
                
            case NonSlipWall:
                *EdgeIndex = 1;
                break;
            case ZeroGrad:
                *EdgeIndex = 1;
                break;
            case Clamped:
                *EdgeIndex = 1;
                break;
            case ClampedDepth:
                *EdgeIndex = 1;
                break;
            case ClampedVel:
                *EdgeIndex = 1;
                break;
            case Flather:
                *EdgeIndex = 1;
                break;
            case NonLinearFlatherDepth:
                *EdgeIndex = 1;
                break;
            case NonLinearFlatherFlow:
                *EdgeIndex = 1;
                break;
            case NonReflectFlux:
                *EdgeIndex = 1;
                break;
            case BottomBoundary:
                *EdgeIndex = 1;
                break;
            case UpperSurfaceBoundary:
                *EdgeIndex = 2;
                break;
        }
}

#endif
