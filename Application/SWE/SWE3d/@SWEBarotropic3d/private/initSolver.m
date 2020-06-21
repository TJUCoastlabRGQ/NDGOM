function [ adv, vis, flux,  PCESolver2d ] = initSolver( physMat )
%INITSOLVER Summary of this function goes here
%   Detailed explanation goes here

integralType = physMat.getOption('integralType');
if (integralType == enumDiscreteIntegral.QuadratureFree)
    adv = NdgQuadFreeStrongFormAdvSolver3d( physMat );
    vis = NdgQuadFreeStrongCentralVisSolver3d( physMat, physMat.varFieldIndex, 1:physMat.Nvar);
    flux = SWEHLLNumFluxSolver3d;
    PCESolver2d = SWEQuadFreeStrongFormPCESolver2d;
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    adv = NdgGaussQuadWeakFormAdvSolver3d( physMat, physMat.meshUnion );
    vis = NdgGaussQuadWeakFormCentralVisSolver3d( physMat, physMat.varFieldIndex, 1:physMat.Nvar);
    flux = SWEHLLNumFluxSolver3d;
    PCESolver2d = SWEGaussQuadWeakFormPCESolver2d( physMat, physMat.mesh2d );
end
% vis = NdgNonVisSolver( physMat );
end



