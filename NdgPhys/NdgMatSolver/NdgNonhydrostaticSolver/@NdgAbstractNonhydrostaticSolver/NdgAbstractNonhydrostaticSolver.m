classdef NdgAbstractNonhydrostaticSolver < handle
    methods
%         function obj = NdgAbstractNonhydrostaticSolver(PhysClass)
%             % Doing nothing
%         end
    end
    
    methods( Access = protected )
        %> @brief function for calculating the Nonhydrostatic pressure
        %> @details
        %> Function to calculate the Nonhydrostatic pressure
        %> @param[in] StiffMatrix the global stiff matrix
        %> @param[in] NonhydrostaticRHS The right hand side of the nonhydrostatic model
        %> @retval[out] NonhydrostaticPressure
        NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
        %> @brief function for assembling the global stiff matrix
        %> @details
        %> Function to calculate the global stiff matrix for the Nonhydrostatic model with conservative correction
        %> @param[in] hm The local water depth
        %> @param[in] hp The adjacent water depth
        %> @param[in] h The water depth
        %> @param[in] physClass The hydrostatic solver
        %> @retval[out] StiffMatrix The global StiffMatrix
        StiffMatrix = matAssembleConservativeGlobalStiffMatrix( obj, hm, hp, h, physClass )
        %> @brief function for calculating the right hand side
        %> @details
        %> Function to calculate the right hand side of the Nonhydrostatic model
        %> @param[in] fphys The fphys field
        %> @param[in] physClass The hydrostatic solver
        %> @retval[out] NonhydrostaticRHS the right hand side of the nonhydrostatic model
        NonhydrostaticRHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass)
        %> @brief function for updating the final flux term and the vertical velocity
        %> @details
        %> Function to calculate the final physics field of the Nonhydrostatic model
        %> @param[in] NonhydrostaticPressure
        %> @param[in] physClass The hydrostatic solver
        %> @retval[out] physClass The final updated hydrostatic solver
        fphys = matUpdateConservativeFinalVelocity(obj, NonhydrostaticPressure, physClass, fphys)
        %> @brief Function to resemble global stiffMatrix and Nonhydrostatic right hand side
        %> @details
        %> Function to get the flux term at the wet dry interface
        %> @param[in] physClass The hydrostatic solver
        %> @param[in] FluxTerm The flux term needs to be corrected
        %> @param[in] mesh The mesh object
        %> @param[out] FluxTerm The flux term corrected
        [StiffMatrix, NonhydrostaticRHS] = matResembleGlobalMatrix(obj, physClass, StiffMatrix, NonhydrostaticRHS);
        %> @brief Function to assemble the wet and dry interface
        %> @details
        %> Function to assemble the wet and dry interface
        %> @param[in] mesh The mesh object
        matAssembleWetDryInterface(obj, mesh)
        %> @brief Function to make the nonhydrostatic correction
        %> @details
        %> Function to make the nonhydrostatic correction
        %> @param[in] physClass The hydrostatic solver
        %> @param[in] fphys The fphys field
        fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys)
        %> @brief Function to set the cell boundary type
        %> @details
        %> Function to  set the cell boundary type to impose the nonhydrostatic related boundary condition
        %> @param[in] mesh The mesh object
        matSetBoundaryType(obj,mesh)
        %> @brief Function to set the characteristic matrix at the initiaizing stage
        %> @details Function to set the characteristic matrix used when assemble the global stiff matrix
        %> @param[in] mesh The mesh object
        matSetInitializeCharacteristicMatrix(obj, PhysClass, mesh)
        %> @brief Function to get the index of the characteristic matrix
        %> @details Function to get the index of the characteristic matrix used to assemble the global stiff sparse matrix
        %> @param[in] PointToCell The matrix PointToCell, used for identifying the cell influenced by points
        %> @param[in] Np the interpolation points of the standard cell
        %> @param[out] JcsStiffMatrix the non-zero element contained in precedent column,  for colume j the number of the non-zero element is JcsGlobalStiffMatrix[j+1] - JcsGlobalStiffMatrix[j]
        %> @param[out] JrsGlobalStiffMatrix the index of the non-zero element in each column, the size of this index is problem dependent, for more information, please refer to mxGetIr
        [JcsStiffMatrix, JrsStiffMatrix] = matGetStiffMatrixIndex(obj, PointToCell, Np)
        %> @brief Function to calculate the characteristic matrix
        %> @details Function to calculate the characteristic matrix in the x direction
        %> @param[in] BoundaryEdge the boundary edge object
        %> @param[in] InnerEdge the inner edge object
        %> @param[in] Variable variable used to calculate the characteristic matrix
        %> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
        %> @param[out] termX the calculated characteristic matrix in x direction
        termX = matCalculateCharacteristicMatrix(obj, BoundaryEdge, InnerEdge, Variable, ftype)
        %> @brief Function to take the local and adjacent face value
        %> @details Function to take the local and adjacent face value
        %> @param[in] fm the local face value
        %> @param[in] fp the adjacent face value
        %> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
        %> @param[out] fm the local face value with wet and dry interface considered
        %> @param[out] fp the adjacent face value with wet and dry interface considered
        [fm, fp] = matGetFaceValue( obj, fm, fp, ftype )
        %> @brief Function to get the adjacent value at the boundary
        %> @details Function to get the adjacent value at the boundary
        %> @param[in] fm the local face value
        %> @param[in] fp the adjacent face value
        %> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
        %> @param[in] EidBoundaryType the boundary index used to impose the zero nonhydrostatic pressure at the open boundary and the zero-grad condition for nonhydrostatic related auxialary variable at the wall
        %> @param[out] fp the adjacent face value at the boundary
        fp = matImposeNonhydroRelatedBoundaryCondition(obj, fm, fp, ftype, EidBoundaryType)
        %> @brief Function to calculate the conservative variable related partial derivative operator
        %> @details Function to calculate the conservative variable related partial derivative operator
        %> @param[in] physClass The hydrostatic solver
        %> @param[in] BoundaryEdge the boundary edge object
        %> @param[in] InnerEdge the inner edge object
        %> @param[in] fphys the physics field
        %> @param[in] ftype the boundary condition at the wet-dry interface
        %> @param[in] index index of the calculated variable
        %> @param[out] termX the discrete operator in x direction
        termX = matCalculateConservativeVariableRelatedMatrixX(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index)
        %> @brief Function to assemble the global sparse stiff matrix
        %> @details Function to assemble the global sparse stiff matrix
        %> @param[in] UpdatedPNPX The updated PNPX
        %> @param[in] UpdatedSPNPX the updated SPNPX
        %> @param[in] UpdatedFNPBX the updated FNPBX
        %> @param[in] fphys the physics field
        %> @param[in] PhysClass the hydrostatic solver
        %> @param[out] StiffMatrix the global sparse stiff matrix
        StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, UpdatedPNPX, UpdatedSPNPX, UpdatedFNPBX, fphys, PhysClass)
        %> @brief Function to assemble the characteristic matrix
        %> @details Function to assemble the characteristic matrix
        %> @param[in] mesh the mesh object
        %> @param[in] index index of the studied point
        %> @param[out] qx the first order derivative of non-hydrostatic prossure with respect to x
        %> @param[out] q2x the second order derivative of non-hydrostatic prossure with respect to x
        %> @param[out] qbx product of non-hydrostatic prossure with bottom gradient in x derection
        %> @param[out] fqbx product of non-hydrostatic prossure with bottom gradient in x derection with flux term function
        [qx, q2x, qbx, fqbx] = matAssembleCharacteristicMatrix(obj, mesh, i)
        %> @brief Function to assemble the point that change status from wet to dry
        %> @details Function to assemble the point that change status from wet to dry
        %> @param[in] Np number of interpolation points
        %> @param[in] Nfp numver of interpolation points per face
        %> @param[in] NewWetDryFace face that change status from wet to wet-dry interface
        %> @param[in] Fmask order of the face interpolation points 
        %> @param[out] TempWetToDryPoint point that change status from wet to dry
        TempWetToDryPoint = matGetTempWetToDryPoint( obj, Np, Nfp, NewWetDryFace, Fmask)
        %> @brief Function to calculate the volume integral in the x direction
        %> @details Function to calculate the volume integral in the x direction
        %> @param[in] mesh the mesh object
        %> @param[in] Variable variable used to calculate the volume integral
        %> @param[out] VolumeIntegralX the volume integral in x direction       
        VolumeIntegralX = matGetVolumeIntegralX(obj, mesh, Variable)
    end
end