classdef NdgNonhydrostaticSolver2d < NdgAbstractNonhydrostaticSolver
    %NDGNONHYDROSTATICSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
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
        function obj = NdgNonhydrostaticSolver2d(PhysClass)
%             obj = obj@NdgAbstractNonhydrostaticSolver(PhysClass);
            mesh = PhysClass.meshUnion(1);
            obj.matSetBoundaryType(mesh);
            obj.NonhydroFmPoint = [];
            obj.NonhydroFpPoint = [];
            [ obj.bx, obj.by ] = obj.matSetBottomGradient(PhysClass.zGrad{1});
            
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
        
    end
    
    methods(  Access = protected, Abstract )
        [ bx, by ] = matSetBottomGradient(obj, zGrad);
        
        matSetInitializeCharacteristicMatrix(obj, physClass, mesh);
        
        [Nq, qx, qy, q2x, q2y] = matAssembleCharacteristicMatrix(obj, mesh, index);
    end
    
    methods( Access = protected )
        %> @brief Function to assemble the global sparse stiff matrix
        %> @details Function to assemble the global sparse stiff matrix
        %> @param[in] UpdatedPNPX The updated PNPX
        %> @param[in] UpdatedSPNPX the updated SPNPX
        %> @param[in] UpdatedFNPBX the updated FNPBX
        %> @param[in] fphys the physics field
        %> @param[in] PhysClass the hydrostatic solver
        %> @param[out] StiffMatrix the global sparse stiff matrix
        StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, UpdatedPNPX, UpdatedPNPY,...
            UpdatedSPNPX, UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY, fphys, PhysClass);
        %> @brief function for calculating the right hand side
        %> @details
        %> Function to calculate the right hand side of the Nonhydrostatic model
        %> @param[in] fphys The fphys field
        %> @param[in] physClass The hydrostatic solver
        %> @retval[out] NonhydrostaticRHS the right hand side of the nonhydrostatic model
        NonhydrostaticRHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass);
        %> @brief function for updating the final flux term and the vertical velocity
        %> @details
        %> Function to calculate the final physics field of the Nonhydrostatic model
        %> @param[in] NonhydrostaticPressure
        %> @param[in] physClass The hydrostatic solver
        %> @retval[out] physClass The final updated hydrostatic solver
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydrostaticPressure, physClass, fphys);
        %> @brief Function to make the nonhydrostatic correction
        %> @details
        %> Function to make the nonhydrostatic correction
        %> @param[in] physClass The hydrostatic solver
        %> @param[in] fphys The fphys field
        fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys);
        %> @brief Function to calculate the characteristic matrix
        %> @details Function to calculate the characteristic matrix in the x direction
        %> @param[in] BoundaryEdge the boundary edge object
        %> @param[in] InnerEdge the inner edge object
        %> @param[in] Variable variable used to calculate the characteristic matrix
        %> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
        %> @param[out] termX the calculated characteristic matrix in x direction        
        [ termX, termY ] = matCalculateCharacteristicMatrix( obj, mesh, BoundaryEdge, InnerEdge, qx, qy, ftype);
        %> @brief Function to calculate the conservative variable related partial derivative operator
        %> @details Function to calculate the conservative variable related partial derivative operator
        %> @param[in] physClass The hydrostatic solver
        %> @param[in] BoundaryEdge the boundary edge object
        %> @param[in] InnerEdge the inner edge object
        %> @param[in] fphys the physics field
        %> @param[in] ftype the boundary condition at the wet-dry interface
        %> @param[in] index index of the calculated variable
        %> @param[out] termX the discrete operator in x direction
        [ termX, termY ] = matCalculateConservativeVariableRelatedMatrix(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index)
        %> @brief Function to calculate the volume integral in the x direction
        %> @details Function to calculate the volume integral in the x direction
        %> @param[in] mesh the mesh object
        %> @param[in] Variable variable used to calculate the volume integral
        %> @param[out] VolumeIntegralX the volume integral in x direction       
        [ VolumeIntegralX, VolumeIntegralY ] = matVolumeIntegral(obj, mesh, VariableX, VariableY )
        
         matAssemblePointToCellInformation(obj, K, Np, PNPX,PNPY, SPNPX, SPNPY, NPBX, NPBY, FNPBX, FNPBY, NP )
         
         [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index)
    end
    
end

