classdef NdgNonhydrostaticSolver2d < NdgAbstractNonhydrostaticSolver
    %NDGNONHYDROSTATICSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %> P stand for partial derivative
        %> NP stands for Non-hydrostatic pressure
        %> X and Y stands for x direction and y direction, respectively
        %> S stands for second order
        %> EidBoundaryType Index used to impose the zero non-hydrostatic boundary conditions at the open boundary,
        %> and to impose the zero grad non-hydrostatic boundary conditions at the wall
        %> boundary conditions at the open boundary, and zero boundary condition at the wall
        %> ZeroFluxBoundary Index used to assemble the wet-dry interface information
        %> ZeroFluxBoundaryIndex Index used to count the wet-dry interface number
        %> AdjacentDryCellAndFace matrix used to record the adjacent dry cell and the corresponding face index
        %> dt time step
        %> JcsGlobalStiffMatrix the non-zero element contained in precedent column,  for colume j the number of the non-zero element is JcsGlobalStiffMatrix[j+1] - JcsGlobalStiffMatrix[j]
        %> JrsGlobalStiffMatrix the index of the non-zero element in each column, the size of this index is problem dependent, for more information, please refer to mxGetIr
        %> HBx the partial derivative of H + 2b with respect to the x direction
        %> HBy the partial derivative of H + 2b with respect to the y direction
        %> H2Bx the second order derivative of H + 2b with respect to the x direction
        %> H2By the second order derivative of H + 2b with respect to the y direction 
        %> fhx the partial derivative of H with respect to the x direction
        %> fhy the partial derivative of H with respect to the y direction
        %> hux the partial derivative of Hu with respect to the x direction
        %> hvy the partial derivative of Hv with respect to the y direction        
        TempPNPX
        PNPX
        TempPNPY
        PNPY
        TempSecondOrderTerm
        SecondOrderTerm
        NP
        WetCellIndex
        TempZeroFluxBoundary
        ZeroFluxBoundary
        ZeroFluxBoundaryIndex
        WetDryFaceOrder
%         NonhydroFmPoint
%         NonhydroFpPoint
        EidBoundaryType
        AdjacentDryCellAndFace
        WetNum
        dt
        WetDryPoint
        TempWetDryPoint
        JcsGlobalStiffMatrix
        JrsGlobalStiffMatrix
        Tau
        edgeType
        HBx
        HBy
        H2Bx
        H2By
        fhx
        fhy
        hux
        hvy
    end
    
    methods
        function obj = NdgNonhydrostaticSolver2d(PhysClass)
%             obj = obj@NdgAbstractNonhydrostaticSolver(PhysClass);
            mesh = PhysClass.meshUnion(1);
            obj.matSetBoundaryType(mesh);
%             obj.NonhydroFmPoint = [];
%             obj.NonhydroFpPoint = [];
            
            obj.matCalculatePenaltyParameter( mesh );
            obj.matAssembleElementBoundaryCondition( mesh );

            obj.matSetInitializeCharacteristicMatrix(PhysClass, mesh);
            obj.matAssemblePointToCellInformation(mesh.K, mesh.cell.Np, obj.PNPX, obj.PNPY,...
                obj.SecondOrderTerm, obj.NP);
            
            obj.ZeroFluxBoundaryIndex = 0;
            obj.ZeroFluxBoundary = ones(0,2);
            obj.TempZeroFluxBoundary = ones(0,2);
            obj.AdjacentDryCellAndFace = [];
            obj.WetDryPoint = [];
            obj.TempWetDryPoint = [];
        end
        
        function fphys = NdgConservativeNonhydrostaticUpdata(obj, PhysClass, fphys, deltatime)
            obj.dt = deltatime;
            fphys = obj.matNdgConservativeNonhydrostaticUpdata( PhysClass, fphys );
        end
        
    end
    
    methods( Access = public, Abstract )
       evaluateNonhydroRHS(obj, PhysClass, fphys);
        
    end
    
    methods(  Access = protected, Abstract )
%         [ bx, by ] = matSetBottomGradient(obj, zGrad);
        
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
            UpdatedSPNPX, UpdatedSPNPY, fphys, PhysClass);
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
        
         matAssemblePointToCellInformation(obj, K, Np, PNPX,PNPY, SPNPX, SPNPY, NP )
        %> @brief Function to calculate the volume integral in the x direction
        %> @details Function to calculate the volume integral in the x direction
        %> @param[in] mesh the mesh object
        %> @param[in] Variable variable used to calculate the volume integral
        %> @param[out] VolumeIntegralX the volume integral in x direction           
         [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index)
         
         [TempHBx, TempHBy] = matCalculateLDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, variable)
         
         [ H2Bx, H2By ] = matCalculateLDGSecondOrderVariable( obj, mesh, BoundaryEdge, InnerEdge, variable, variablex, variabley )
                  
         matCalculateFphysDerivative(obj, mesh, fphys, physClass)
    end
    
end

