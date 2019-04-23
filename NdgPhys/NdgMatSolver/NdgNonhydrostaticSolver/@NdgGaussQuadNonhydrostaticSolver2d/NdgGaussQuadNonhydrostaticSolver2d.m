classdef NdgGaussQuadNonhydrostaticSolver2d < NdgNonhydrostaticSolver2d & ...
        NdgGaussQuadWeakFormSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver2d( physClass, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( physClass, meshUnion );
        end
    end
    
    methods( Access = protected )
        [qx, qy, q2x, q2y, qbx, qby, fqbx, fqby, Nonhydrop] = matAssembleCharacteristicMatrix(obj, mesh, index);
        
        [ tempqx, tempqy ]  = matCalculateCharacteristicMatrix( obj, mesh, BoundaryEdge, InnerEdge, VariableX, VariableY, ftype);
        
        [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, PhysClass, BoundaryEdge, InnerEdge, fphys, ftype, index);        
        
        RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass);
                        
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys);
        
        [ VolumeIntegralX, VolumeIntegralY ] = matVolumeIntegral(obj, mesh, VariableX, VariableY);
        
    end
    
end

