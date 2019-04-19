function [ adv, vis, wind, flux, Friction, PCESolver2d ] = initSolver( physMat )
%INITSOLVER Summary of this function goes here
%   Detailed explanation goes here

integralType = physMat.getOption('integralType');
if (integralType == enumDiscreteIntegral.QuadratureFree)
    adv = NdgQuadFreeStrongFormAdvSolver3d( physMat );
    vis = NdgQuadFreeStrongCentralVisSolver3d( physMat, physMat.varFieldIndex, 1:physMat.Nvar);
    wind = NdgQuadFreeWindSolver3d();
    flux = SWELFNumFluxSolver3d;
    Friction = NdgQuadFreeFrictionSolver3d;
    PCESolver2d = SWEQuadFreeStrongFormPCESolver2d;
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    adv = NdgGaussQuadWeakFormAdvSolver3d( physMat, physMat.meshUnion );
    vis = NdgGaussQuadWeakFormCentralVisSolver3d( physMat, physMat.varFieldIndex, 1:physMat.Nvar);
    wind = NdgGaussQuadWindSolver3d( physMat, physMat.meshUnion );
    flux = SWELFNumFluxSolver3d;
    Friction = NdgGaussQuadFrictionSolver3d( physMat, physMat.meshUnion );    
    PCESolver2d = SWEGaussQuadWeakFormPCESolver2d( physMat, physMat.mesh2d );
end
% vis = NdgNonVisSolver( physMat );
end



