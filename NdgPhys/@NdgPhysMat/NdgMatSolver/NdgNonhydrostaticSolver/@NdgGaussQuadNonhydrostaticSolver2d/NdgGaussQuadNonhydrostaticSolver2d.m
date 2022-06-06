classdef NdgGaussQuadNonhydrostaticSolver2d < NdgNonhydrostaticSolver2d & ...
        NdgGaussQuadWeakFormSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver2d( physClass, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( physClass, meshUnion );
            obj = obj@NdgNonhydrostaticSolver2d(physClass);
        end
        
        function evaluateNonhydroRHS(obj, PhysClass, fphys)
%             mesh = PhysClass.meshUnion(1);
%             BoundaryEdge = mesh.BoundaryEdge;
%             InnerEdge = mesh.InnerEdge;
%             [ HBx, HBy ] = obj.matCalculateCharacteristicMatrix( mesh,  BoundaryEdge, ...
%                 InnerEdge, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]),  num2cell(fphys{1}(:,:,1) +...
%                 2 * fphys{1}(:,:,4),[1 2]), enumNonhydroBoundaryCondition.Zero);
%             
%             [ px, py ] = obj.matCalculateCharacteristicMatrix( mesh,  BoundaryEdge, ...
%                 InnerEdge, num2cell(fphys{1}(:,:,7),[1 2]),  num2cell(fphys{1}(:,:,7),[1 2]),...
%                 enumNonhydroBoundaryCondition.Zero);
%             
%             PhysClass.frhs{1}(:,:,2) = PhysClass.frhs{1}(:,:,2) - fphys{1}(:,:,1) .* px - fphys{1}(:,:,7) .* HBx;
%             PhysClass.frhs{1}(:,:,3) = PhysClass.frhs{1}(:,:,3) - fphys{1}(:,:,1) .* py - fphys{1}(:,:,7) .* HBy;
%             PhysClass.frhs{1}(:,:,4) = PhysClass.frhs{1}(:,:,4) + 2 * fphys{1}(:,:,7);
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

         matReconstructStiffmatrixRelatedMatrix( obj, physClass)
        
    end
    
end

