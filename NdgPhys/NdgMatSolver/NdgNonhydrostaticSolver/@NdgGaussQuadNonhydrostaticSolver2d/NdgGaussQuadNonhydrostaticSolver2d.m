classdef NdgGaussQuadNonhydrostaticSolver2d < NdgNonhydrostaticSolver2d & ...
        NdgGaussQuadWeakFormSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver2d( physClass, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( physClass, meshUnion );
            obj = obj@NdgNonhydrostaticSolver2d(physClass);
        end
    end
    
    methods( Access = protected )
                
        fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys);
        
        matSetInitializeCharacteristicMatrix(obj, physClass, mesh);
        
        [Nq, qx, qy, q2x, q2y] = matAssembleCharacteristicMatrix(obj, mesh, index);
        
        [ tempqx, tempqy ]  = matCalculateCharacteristicMatrix( obj, mesh, BoundaryEdge, InnerEdge, VariableX, VariableY, ftype);
        
        [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, physClass, mesh, fphys, ftype, index);        
        
        RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass);
                        
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys);
        
        [ VolumeIntegral ] = matVolumeIntegral( obj, Variable );
        
        [ VolumeIntegralX, VolumeIntegralY ] = matFluxVolumeIntegral(obj, VariableX, VariableY);
        
        [ VolumeIntegralX, VolumeIntegralY ] = matEvaluateConservativeVariableTotalIntegral( obj, physClass, mesh, ftype, Variable, VolumeIntegralX, VolumeIntegralY, Index );
        
        [ VolumeIntegralX, VolumeIntegralY ] = matEvaluateNonhydroVariableTotalIntegral( obj, mesh, ftype, VariableX, VariableY, VolumeIntegralX, VolumeIntegralY );
        
        [ bx, by ] = matSetBottomGradient( obj, zGrad);
        
        StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix( obj, PhysClass, UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, UpdatedNP, fphys)

        [UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, UpdatedNP ] = matReconstructStiffmatrixRelatedMatrix( obj, physClass)
        
    end
    
end

