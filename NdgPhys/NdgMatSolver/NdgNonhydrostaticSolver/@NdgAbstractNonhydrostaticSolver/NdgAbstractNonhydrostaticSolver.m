classdef NdgAbstractNonhydrostaticSolver < handle
    methods(Abstract, Access = protected)
        %> @brief Function to assemble the global sparse stiff matrix
        %> @details Function to assemble the global sparse stiff matrix
        %> @param[in] UpdatedPNPX The updated PNPX
        %> @param[in] UpdatedSPNPX the updated SPNPX
        %> @param[in] UpdatedFNPBX the updated FNPBX
        %> @param[in] fphys the physics field
        %> @param[in] PhysClass the hydrostatic solver
        %> @param[out] StiffMatrix the global sparse stiff matrix
        StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, UpdatedPNPX, UpdatedSPNPX, UpdatedFNPBX, fphys, PhysClass)
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
        %> @brief Function to make the nonhydrostatic correction
        %> @details
        %> Function to make the nonhydrostatic correction
        %> @param[in] physClass The hydrostatic solver
        %> @param[in] fphys The fphys field
        fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys)
        %> @brief Function to set the characteristic matrix at the initiaizing stage
        %> @details Function to set the characteristic matrix used when assemble the global stiff matrix
        %> @param[in] mesh The mesh object
        matSetInitializeCharacteristicMatrix(obj, PhysClass, mesh)
        %> @brief Function to calculate the characteristic matrix
        %> @details Function to calculate the characteristic matrix in the x direction
        %> @param[in] BoundaryEdge the boundary edge object
        %> @param[in] InnerEdge the inner edge object
        %> @param[in] Variable variable used to calculate the characteristic matrix
        %> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
        %> @param[out] termX the calculated characteristic matrix in x direction
        matCalculateCharacteristicMatrix ( obj )
        %> @brief Function to assemble the characteristic matrix
        %> @details Function to assemble the characteristic matrix
        %> @param[in] mesh the mesh object
        %> @param[in] index index of the studied point
        %> @param[out] qx the first order derivative of non-hydrostatic prossure with respect to x
        %> @param[out] q2x the second order derivative of non-hydrostatic prossure with respect to x
        %> @param[out] qbx product of non-hydrostatic prossure with bottom gradient in x derection
        %> @param[out] fqbx product of non-hydrostatic prossure with bottom gradient in x derection with flux term function
        matAssembleCharacteristicMatrix( obj );  
        
        matVolumeIntegral(obj)
    end
    
    methods( Access = protected )
        %> @brief function for calculating the Nonhydrostatic pressure
        %> @details
        %> Function to calculate the Nonhydrostatic pressure
        %> @param[in] StiffMatrix the global stiff matrix
        %> @param[in] NonhydrostaticRHS The right hand side of the nonhydrostatic model
        %> @retval[out] NonhydrostaticPressure
        NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
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
        %> @brief Function to set the cell boundary type
        %> @details
        %> Function to  set the cell boundary type to impose the nonhydrostatic related boundary condition
        %> @param[in] mesh The mesh object
        matSetBoundaryType(obj,mesh)
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
        %> @brief Function to assemble the point that change status from wet to dry
        %> @details Function to assemble the point that change status from wet to dry
        %> @param[in] Np number of interpolation points
        %> @param[in] Nfp numver of interpolation points per face
        %> @param[in] NewWetDryFace face that change status from wet to wet-dry interface
        %> @param[in] Fmask order of the face interpolation points 
        %> @param[out] TempWetToDryPoint point that change status from wet to dry
        TempWetToDryPoint = matGetTempWetToDryPoint( obj, Np, Nfp, NewWetDryFace, Fmask)
        %> @brief Function to calculate the numerical flux of the inner edge corresponding to the primitive variable
        %> with LDG method adopted
        %> @details Function to calculate the numerical flux of the inner edge corresponding to the primitive variable
        %> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided 
        %> by Peraire and Perssony (2008)
        %> @param[in] WetDryFaceOrder order of wet-dry interface
        %> @param[in] fluxS the numerical flux without wet-dry interface taking into consideration
        %> @param[in] Nfp number of interpolation points per face
        %> @param[out] fluxS the numerical flux with wet-dry interface taking into consideration      
        fluxS = matGetPrimitiveVariableInnerEdgeFlux( obj, WetDryFaceOrder, fluxS, Nfp)
        %> @brief Function to calculate the numerical flux of the boundary edge corresponding to the primitive variable
        %> with LDG method adopted
        %> @details Function to calculate the numerical flux of the boundary edge corresponding to the primitive variable
        %> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided 
        %> by Peraire and Perssony (2008)
        %> @param[in] vector the directional vector in x or y direction
        %> @param[in] Um value of the local non-hydrostatic pressure
        %> @param[out] fluxS the numerical flux in direction x or y 
        fluxS = matGetPrimitiveVariableBoundaryEdgeFlux( obj, vector, Um) 
        %> @brief Function to calculate the numerical flux of the inner edge corresponding to the auxialary variable
        %> with LDG method adopted
        %> @details Function to calculate the numerical flux of the inner edge corresponding to the auxialary variable
        %> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided 
        %> by Peraire and Perssony (2008)
        %> @param[in] fluxS the numerical flux without wet-dry interface taking into consideration
        %> @param[in] Um value of the local non-hydrostatic pressure
        %> @param[in] Up value of the adjacent non-hydrostatic pressure
        %> @param[in] Sigmam value of the local auxialary pressure
        %> @param[in] Sigmap value of the adjacent auxialary pressure        
        %> @param[in] vector the directional vector in x or y direction
        %> @param[out] fluxS the numerical flux with wet-dry interface taking into consideration           
        fluxS = matGetAuxialaryVariableInnerEdgeFlux( obj, fluxS, Um, Up, Sigmam, Sigmap, vector)
        %> @brief Function to calculate the numerical flux of the boundary edge corresponding to the auxialary variable
        %> with LDG method adopted
        %> @details Function to calculate the numerical flux of the boundary edge corresponding to the auxialary variable
        %> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided 
        %> by Peraire and Perssony (2008)
        %> @param[in] Um value of the local non-hydrostatic pressure
        %> @param[in] Sigmam value of the local auxialary pressure
        %> @param[in] vector the directional vector in x or y direction
        %> @param[out] fluxS the numerical flux with boundary contion taking into consideration           
        fluxS = matGetAuxialaryVariableBoundaryEdgeFlux( obj,  Um, Sigmam, vector)
        %> @brief Function to get the variable value and the numerical flux at the wet-dry interface
        %> @details Function to get the variable value and the numerical flux at the wet-dry interface
        %> @param[in] Um value of the local cell
        %> @param[in] Up value of the adjacent cell
        %> @param[in] fluxS the flux without wet-dry interface considered
        %> @param[out] Um value of the local cell with wet-dry interface considered   
        %> @param[out] Up value of the adjacent cell with wet-dry interface considered
        %> @param[in] fluxS the flux with wet-dry interface considered
        [Um, Up, fluxS] = matGetWetDryFaceVariableAndFlux( obj, Um, Up, fluxS)
    end
end