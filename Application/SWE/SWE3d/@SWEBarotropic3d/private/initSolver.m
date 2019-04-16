function [ adv, vis, wind, flux, Friction, PCESolver2d ] = initSolver( physMat )
%INITSOLVER Summary of this function goes here
%   Detailed explanation goes here

integralType = physMat.getOption('integralType');
if (integralType == enumDiscreteIntegral.QuadratureFree)
    adv = NdgQuadFreeStrongFormAdvSolver3d( physMat );
    vis = NdgQuadFreeStrongCentralVisSolver3d( physMat, physMat.varFieldIndex, 1:numel(physMat.varFieldIndex));
    wind = NdgQuadFreeWindSolver3d();
    flux = SWELFNumFluxSolver3d;
    Friction = NdgQuadFreeFrictionSolver3d;
    PCESolver2d = SWEQuadFreeStrongFormPCESolver2d;
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    adv = NdgGaussQuadWeakFormAdvSolver3d( physMat );
    vis = NdgGaussQuadWeakFormVisSolver3d( physMat, physMat.Nvar, 1:numel(physMat.Nvar));
    wind = NdgGaussQuadWindSolver3d();
    flux = SWELFNumFluxSolver3d;
    Friction = NdgGaussQuadFrictionSolver3d;    
    PCESolver2d = SWEGaussQuadPCESolver2d;
end
% vis = NdgNonVisSolver( physMat );
end



