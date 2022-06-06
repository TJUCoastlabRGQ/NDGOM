classdef NdgAbstractNonhydrostaticSolver < handle
    
    properties
        
        EidBoundaryType
        
        Tau
        
        edgeType
        
        PNPX
        
        TempPNPX
        
        SecondOrderTerm
        
        TempSecondOrderTerm
        
        NP
        
        JcsGlobalStiffMatrix
        
        JrsGlobalStiffMatrix
        
    end
    
    properties
        
        dt
        
        WetDryCell
        
        HBx
        
        H2Bx
        
        fhx
        
        hux
        
    end
    
    methods
        function obj = NdgAbstractNonhydrostaticSolver(PhysClass)
            mesh = PhysClass.meshUnion(1);
            obj.matSetBoundaryType(mesh);
            obj.matCalculatePenaltyParameter( mesh );
            obj.matAssembleElementBoundaryCondition( mesh );
        end
        
        function fphys = NdgConservativeNonhydrostaticUpdata(obj, PhysClass, fphys, deltatime)
            obj.dt = deltatime;
            fphys = obj.matNdgConservativeNonhydrostaticUpdata( PhysClass, fphys );   
        end
        
        function matClearGlobalMemory(obj)
            %doing nothing
        end
        
    end
    
    methods(Abstract, Access = protected)
        %> @brief Function to set up the matrix corresponding the Laplace operator and
        %> the first order operator related to nonhydrostatic pressure
        %> @details
        %> Function to set up the matrix corresponding the Laplace operator
        %> and the first order operator related to nonhydrostatic pressure
        %> @param[in] physClass the physical object set up
        %> @param[in] mesh The mesh object
        matSetInitializeCharacteristicMatrix(obj, physClass, mesh)
        %> @brief Function to make the nonhydrostatic correction
        %> @details
        %> Function to make the nonhydrostatic correction
        %> @param[in] physClass The physical object set up
        %> @param[in] fphys The fphys field
        fphys = matNdgConservativeNonhydrostaticUpdata( obj, physClass, fphys )
        %> @brief Function to reassemble the characteristic matrix
        %> @details
        %> Function to reassemble the characteristic matrix with wetting and drying considered
        %> @param[in] mesh The mesh object
        matReconstructStiffmatrixRelatedMatrix(obj , mesh)
        %> @brief Function to calculate the physical variable related partial derivative
        %> @details
        %> Function to calculate the physical variable related partial derivative
        %> @param[in] mesh The mesh object
        %> @param[in] fphys The fphys field     
        %> @param[in] physClass The physical object set up   
        matCalculateFphysDerivative( obj, mesh, fphys, physClass )
        %> @brief Function to assemble the global sparse stiff matrix
        %> @details Function to assemble the global sparse stiff matrix
        %> @param[in] fphys the physics field
        %> @param[out] StiffMatrix the global sparse stiff matrix        
        StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix( obj, fphys )
        %> @brief Function to calculate the right hand side of the non-hydrostatic pressure
        %> @details
        %> Function to calculate the right hand side of the non-hydrostatic pressure
        %> @param[in] fphys The fphys field     
        %> @param[in] physClass The physical object set up
        %> @param[out] NonhydrostaticRHS the right hand side corresponding to the non-hydrostatic part
        NonhydrostaticRHS  = matEvaluateConservativeNonhydrostaticRHS( obj, fphys, physClass )
        %> @brief function for obtaining the divergence free velocity
        %> @details
        %> Function to calculate the divergence-free velocity field
        %> @param[in] NonhydrostaticPressure the non-hydrostatic pressure
        %> @param[in] physClass The physical object set up
        %> @param[in] fphys The fphys field 
        %> @retval[out] fphys The divergence-free fphys field        
        fphys = matUpdateConservativeFinalVelocity( obj, NonhydrostaticPressure , physClass, fphys)        
    end
    
    methods( Access = protected )
        %> @brief Function to set the cell boundary type
        %> @details
        %> Function to  set the cell boundary type to impose the nonhydrostatic related boundary condition
        %> @param[in] mesh The mesh object
        matSetBoundaryType(obj, mesh)
        %> @brief Function to calculate the penalty parameter for the IPDG method
        %> @details
        %> Function to calculate the penalty parameter for the IPDG method
        %> @param[in] mesh The mesh object
        matCalculatePenaltyParameter( obj, mesh )
        %> @brief Function to assemble the boundary condition of each edge for a give mesh cell
        %> @details
        %> Function to assemble the boundary condition of each edge for a give mesh cell
        %> @param[in] mesh The mesh object
        matAssembleElementBoundaryCondition( obj, mesh )
        %> @brief Function to aspect the influence scope of a given study interpolation point
        %> @details
        %> Function to aspect the influence scope of a given study interpolation point,
        %> this could help to know advance the structure of the final stiff matrix, and
        %> helps to adopt openmp for accelaration
        %> @param[in] K The cell number of the studied mesh
        %> @param[in] Np The number of interpolation points per cell
        matAssemblePointToCellInformation( obj, K, Np)
        %> @brief Function to resemble global stiffMatrix and Nonhydrostatic right hand side
        %> @details
        %> Function to get the flux term at the wet dry interface
        %> @param[in] physClass The physical object set up   
        %> @param[in] FluxTerm The flux term needs to be corrected
        %> @param[in] mesh The mesh object
        %> @param[out] FluxTerm The flux term corrected        
        [StiffMatrix, NonhydrostaticRHS] = matResembleGlobalMatrix(obj, mesh, StiffMatrix, NonhydrostaticRHS)
        %> @brief function for calculating the Nonhydrostatic pressure
        %> @details
        %> Function to calculate the Nonhydrostatic pressure
        %> @param[in] StiffMatrix the global stiff matrix
        %> @param[in] NonhydrostaticRHS The right hand side of the nonhydrostatic model
        %> @retval[out] NonhydrostaticPressure
        NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
        %> @brief Function to get the numerical flux of the inner edge corresponding to the primitive variable
        %> with LDG wet-dry condition considered
        %> @details Function to get the numerical flux of the inner edge corresponding to the primitive variable
        %> with LDG wet-dry condition considered
        %> @param[in] WetDryFaceOrder order of wet-dry interface
        %> @param[in] fluxS the numerical flux without wet-dry interface taking into consideration
        %> @param[in] Nfp number of interpolation points per face
        %> @param[out] fluxS the numerical flux with wet-dry interface taking into consideration        
        fluxS = matGetPrimitiveVariableInnerEdgeFlux(  obj, WetDryFaceOrder, fluxS, Nfp)
        %> @brief Function to calculate the numerical flux of the boundary edge corresponding to the primitive variable
        %> @details Function to calculate the numerical flux of the boundary edge corresponding to the primitive variable
        %> @param[in] vector the directional vector in x or y direction
        %> @param[in] Um value of the local non-hydrostatic pressure
        %> @param[out] fluxS the numerical flux in direction x or y
        fluxS = matGetPrimitiveVariableBoundaryEdgeFlux( obj, vector, Um)
    end
end