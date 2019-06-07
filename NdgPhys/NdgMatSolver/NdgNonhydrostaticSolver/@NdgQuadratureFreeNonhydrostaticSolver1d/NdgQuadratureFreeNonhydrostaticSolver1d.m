classdef NdgQuadratureFreeNonhydrostaticSolver1d < NdgNonhydrostaticSolver1d
    %NDGQUADRATUREFREENONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        
        function obj = NdgQuadratureFreeNonhydrostaticSolver1d(PhysClass)
            
            obj = obj@NdgNonhydrostaticSolver1d(PhysClass);
            
        end
        
    end
    
    methods(Access=protected)
        
        [qx, q2x, Np] = ...
            matAssembleCharacteristicMatrix(obj, mesh, index);
        
        RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass);
        
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys);
        
        [ termX ] = matCalculateCharacteristicMatrix(obj, mesh,  BoundaryEdge, InnerEdge, VariableX, ftype)
        
        [ VolumeIntegralX ] = matVolumeIntegral(obj, mesh, VariableX );
        
        [ termX ] = matCalculateConservativeVariableRHSMatrix( obj, PhysClass, BoundaryEdge, InnerEdge, fphys, ftype, index);
        
        matSetInitializeCharacteristicMatrix(obj, physClass, mesh);
        
        [ qx ] = matCalculateLDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable);
        
        [ q2x ] = matCalculateLDGSecondOrderVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable, VariableX );
        
        matCalculateLDGPenaltyParameter(obj, mesh);
        
        matCalculateFphysDerivative(obj, mesh, fphys, physClass);
        
        [ UpwindedTermX ] = matCalculateUpwindedFphysDerivative( obj, mesh, fphys, physClass, variableX );
        
        [ DownwindedTermX ] = matCalculateDownwindedFphysDerivative( obj, mesh, fphys, physClass, variableX );
        
    end    
    
end

