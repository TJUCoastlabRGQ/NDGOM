classdef NdgQuadratureFreeNonhydrostaticSolver2d < NdgNonhydrostaticSolver2d
    
    %NDGNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        
    end
    
    methods
        
        function obj = NdgQuadratureFreeNonhydrostaticSolver2d(PhysClass)
            obj = obj@NdgNonhydrostaticSolver2d(PhysClass);
        end
        
        function evaluateNonhydroRHS(obj, PhysClass, fphys)
            mesh = PhysClass.meshUnion(1);
            NonhydroVolumeflux = 1/2 * fphys{1}(:,:,7) .* fphys{1}(:,:,1);
            [ NqHx , NqHy ] = obj.matCalculateCharacteristicMatrix( mesh,  mesh.BoundaryEdge, mesh.InnerEdge,...
                num2cell(NonhydroVolumeflux,[1 2]), num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);
%             NqHx = mesh.rx .* ( mesh.cell.Dr * NonhydroVolumeflux ) + mesh.sx .* ( mesh.cell.Ds * NonhydroVolumeflux);
%             NqHy = mesh.ry .* ( mesh.cell.Dr * NonhydroVolumeflux ) + mesh.sy .* ( mesh.cell.Ds * NonhydroVolumeflux);
            PhysClass.frhs{1}(:,:,2) = PhysClass.frhs{1}(:,:,2) - NqHx - fphys{1}(:,:,7) .* obj.bx;
            PhysClass.frhs{1}(:,:,3) = PhysClass.frhs{1}(:,:,3) - NqHy - fphys{1}(:,:,7) .* obj.by;
            PhysClass.frhs{1}(:,:,4) = PhysClass.frhs{1}(:,:,4) + fphys{1}(:,:,7);
        end
        
    end
    
    methods(Access=protected)
        
        [qx, qy, q2x, q2y, Np] = ...
            matAssembleCharacteristicMatrix(obj, mesh, index);
        
        RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass);
        
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys);
        
        [ termX, termY] = matCalculateCharacteristicMatrix(obj, mesh,  BoundaryEdge, InnerEdge, VariableX, VariableY, ftype)
        
        [ VolumeIntegralX, VolumeIntegralY ] = matVolumeIntegral(obj, mesh, VariableX, VariableY);
        
        [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, PhysClass, BoundaryEdge, InnerEdge, fphys, ftype, index);
                
        matSetInitializeCharacteristicMatrix(obj, physClass, mesh);
        
        [qx, qy] = matCalculateLDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable);
        
        [q2x, q2y] = matCalculateLDGSecondOrderVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable, VariableX, VariableY );
        
    end
    methods
        %> Functions following are used for testing purpose
        function getWetDryInterface(obj, mesh)
            obj.matAssembleWetDryInterface(mesh);
        end
        function Outer = getFaceOuterValue(obj,mesh, Inner, Outer)
            Outer = obj.matGetFaceOuterValue( mesh, Inner, Outer);
        end
        function fluxterm = GetFluxTerm(obj, mesh, fluxterm)
            fluxterm = obj.matGetFluxTerm(mesh, fluxterm);
        end
        function [fm, fp] = GetFaceValue( obj, fm, fp, ftype)
            [fm, fp] = obj.matGetFaceValue(fm, fp, ftype);
        end
        
        function [qx, qy] = GetCharacteristicMatrix( obj, mesh, gmatx, gmaty, ftype)
            BoundaryEdge = mesh.BoundaryEdge; InnerEdge = mesh.InnerEdge;
            qx =  obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, gmatx, ftype);
            qy =  obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, gmaty, ftype);
        end
        
        function fp = GetBoundaryValue(obj, fm, fp, ftype)
            fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, ftype, obj.EidBoundaryType);
        end
        
        function StiffMatrix = GetGlobalStiffMatrix(obj, PNPX, PNPY, SPNPX, SPNPY, FNPBX, FNPBY, fphys, PhysClass)
            StiffMatrix = obj.matAssembleConservativeGlobalSparseStiffMatrix( PNPX, PNPY,...
                SPNPX, SPNPY, FNPBX, FNPBY, fphys, PhysClass);
        end
    end
end