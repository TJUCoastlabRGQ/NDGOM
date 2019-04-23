classdef NdgQuadratureFreeNonhydrostaticSolver2d < NdgNonhydrostaticSolver2d

    %NDGNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties

    end
    
    methods(Access=protected)
                        
        [qx, qy, q2x, q2y, qbx, qby, fqbx, fqby, Np] = ...
            matAssembleCharacteristicMatrix(obj, mesh, i);
                
        StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, UpdatedPNPX, UpdatedPNPY,...
            UpdatedSPNPX, UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY, fphys, PhysClass);
        
        RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass);
        
        termY = matCalculateConservativeVariableRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index);
                
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys);
        
        [ termX, termY] = matCalculateCharacteristicMatrix(obj, mesh,  BoundaryEdge, InnerEdge, VariableX, VariableY, ftype)
        
        [ VolumeIntegralX, VolumeIntegralY ] = matVolumeIntegral(obj, mesh, VariableX, VariableY);
        
        [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, PhysClass, BoundaryEdge, InnerEdge, fphys, ftype, index);
        
    end
    
    methods
        
        function obj = NdgQuadratureFreeNonhydrostaticSolver2d(PhysClass)
            obj = obj@NdgNonhydrostaticSolver2d(PhysClass);
        end        
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
        
        function matrixY = GetCharacteristicMatrixY(obj, mesh, BoundaryEdge, InnerEdge, Variable, ftype)
            matrixY = obj.matCalculateCharacteristicMatrixY(mesh, BoundaryEdge, InnerEdge, Variable, ftype);
        end 
        
        function matrixX = GetConservativeVariableRelatedMatrixX(obj, BoundaryEdge, InnerEdge, fphys, ftype, index)
            matrixX = obj.matCalculateConservativeVariableRelatedMatrixX(obj.phys, BoundaryEdge, InnerEdge, fphys, ftype, index);
        end
        
        function matrixY = GetConservativeVariableRelatedMatrixY(obj, BoundaryEdge, InnerEdge, fphys, ftype, index)
            matrixY = obj.matCalculateConservativeVariableRelatedMatrixY(obj.phys, BoundaryEdge, InnerEdge, fphys, ftype, index);
        end 
        
        function StiffMatrix = GetGlobalStiffMatrix(obj, PNPX, PNPY, SPNPX, SPNPY, FNPBX, FNPBY, fphys, PhysClass)
                 StiffMatrix = obj.matAssembleConservativeGlobalSparseStiffMatrix( PNPX, PNPY,...
            SPNPX, SPNPY, FNPBX, FNPBY, fphys, PhysClass);
        end        
    end
end