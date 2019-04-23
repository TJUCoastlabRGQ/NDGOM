classdef NdgQuadratureFreeNonhydrostaticSolver2d < NdgQuadratureFreeNonhydrostaticSolver

    %NDGNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %> NonhydroFlag used to indicate whether nonhydrostatic correction or not
        %> P stand for partial derivative
        %> NP stands for Non-hydrostatic pressure
        %> X and Y stands for x direction and y direction, respectively
        %> S stands for second order
        %> BX and BY stand for topography gradient in x and y direction, respectively
        %> F stands for flux
        %> EidBoundaryType Index used to impose the zero non-hydrostatic boundary conditions at the open boundary,
        %> and to impose the zero grad non-hydrostatic boundary conditions at the wall
        %> boundary conditions at the open boundary, and zero boundary condition at the wall
        %> ZeroFluxBoundary Index used to assemble the wet-dry interface information
        %> ZeroFluxBoundaryIndex Index used to count the wet-dry interface number
        %> AdjacentDryCellAndFace matrix used to record the adjacent dry cell and the corresponding face index
        %> bx the bottom topography gradient in the x direction, with the boundary gradient term averaged
        %> by the bottom topography gradient in the y direction, with the boundary gradient term averaged
        %> dt time step
        %> JcsGlobalStiffMatrix the non-zero element contained in precedent column,  for colume j the number of the non-zero element is JcsGlobalStiffMatrix[j+1] - JcsGlobalStiffMatrix[j]
        %> JrsGlobalStiffMatrix the index of the non-zero element in each column, the size of this index is problem dependent, for more information, please refer to mxGetIr
        NonhydroFlag
        TempPNPX
        PNPX
        TempPNPY
        PNPY
        TempSPNPX
        SPNPX
        TempSPNPY
        SPNPY
        NPBX
        NPBY
        NP
        TempFNPBX
        FNPBX
        TempFNPBY
        FNPBY
        WetCellIndex
        TempZeroFluxBoundary
        ZeroFluxBoundary
        ZeroFluxBoundaryIndex
        NonhydroFmPoint
        NonhydroFpPoint
        EidBoundaryType
        AdjacentDryCellAndFace
        WetNum
        bx
        by
        dt
        WetDryPoint
        TempWetDryPoint
        JcsGlobalStiffMatrix
        JrsGlobalStiffMatrix
    end
    
    methods
        function obj = NdgQuadratureFreeNonhydrostaticSolver2d(PhysClass)
%             obj = obj@NdgAbstractNonhydrostaticSolver(PhysClass);
            mesh = PhysClass.meshUnion(1);
            obj.matSetBoundaryType(mesh);
            obj.NonhydroFmPoint = [];
            obj.NonhydroFpPoint = [];
            obj.bx = PhysClass.zGrad{1}(:,:,1);
            obj.by = PhysClass.zGrad{1}(:,:,2);
            obj.matSetInitializeCharacteristicMatrix(PhysClass, mesh);
            obj.matAssemblePointToCellInformation(mesh.K, mesh.cell.Np, obj.PNPX, obj.PNPY, obj.SPNPX, obj.SPNPY,...
                obj.NPBX,obj.NPBY,obj.FNPBX, obj.FNPBY,obj.NP);
            obj.ZeroFluxBoundaryIndex = 0;
            obj.ZeroFluxBoundary = ones(0,2);
            obj.TempZeroFluxBoundary = ones(0,2);
            obj.AdjacentDryCellAndFace = [];
            obj.WetDryPoint = [];
            obj.TempWetDryPoint = [];
            %             obj.NonhydroFlag = 0;
        end
        
        function fphys = NdgConservativeNonhydrostaticUpdata(obj, PhysClass, fphys, deltatime)
            obj.dt = deltatime;
            fphys = obj.matNdgConservativeNonhydrostaticUpdata( PhysClass, fphys );
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
    
    methods(Access=protected)
        
        fphys = matNdgConservativeNonhydrostaticUpdata( obj, PhysClass, fphys );
        
        matSetInitializeCharacteristicMatrix(obj, PhysClass, mesh);
        
        [qx, qy, q2x, q2y, qbx, qby, fqbx, fqby, Np] = ...
            matAssembleCharacteristicMatrix(obj, mesh, i);
        
        matAverageBoundaryTopoGrad(obj, PhysClass, mesh);
        
        StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, UpdatedPNPX, UpdatedPNPY,...
            UpdatedSPNPX, UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY, fphys, PhysClass);
        
        RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass);
        
        termY = matCalculateConservativeVariableRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index);
        
        termY = matCalculateCharacteristicMatrixY(obj, mesh, BoundaryEdge, InnerEdge, Variable, ftype);
        
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys);
        
       [ LDGqx, LDGqy ] = matAlternatingUpwindLDGFlux( obj, mesh, BoundaryEdge, InnerEdge, gmat, ftype); 
       [ LDGq2x, LDGq2y] = matAlternatingDownwindLDGFlux( obj, mesh, BoundaryEdge, InnerEdge, LDGqx, ftype);
    end
end