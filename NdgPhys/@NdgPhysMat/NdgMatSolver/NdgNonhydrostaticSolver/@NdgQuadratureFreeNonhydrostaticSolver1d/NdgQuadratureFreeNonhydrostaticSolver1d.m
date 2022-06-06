classdef NdgQuadratureFreeNonhydrostaticSolver1d < NdgNonhydrostaticSolver1d
    %NDGQUADRATUREFREENONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    %NDGNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        
    end
    
    methods
        
        function obj = NdgQuadratureFreeNonhydrostaticSolver1d(PhysClass)
            obj = obj@NdgNonhydrostaticSolver1d(PhysClass);
            obj.matAssemblePointToCellInformation( PhysClass.meshUnion(1).K, ...
                PhysClass.meshUnion(1).cell.Np );
        end
        
        function evaluateNonhydroRHS(obj, physClass, fphys)
%             mesh = physClass.meshUnion(1);
%             BoundaryEdge = mesh.BoundaryEdge;
%             InnerEdge = mesh.InnerEdge;
%             [ HBx ] = obj.matCalculateUpwindedFphysDerivative( mesh,  fphys, ...
%                 physClass, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,3),[1 2]));
%             [ px ] = obj.matCalculateIPDGAuxialaryVariable( mesh,  BoundaryEdge, ...
%                 InnerEdge, num2cell(fphys{1}(:,:,6),[1 2]));
%             physClass.frhs{1}(:,:,2) = physClass.frhs{1}(:,:,2) - fphys{1}(:,:,1) .* px - fphys{1}(:,:,6) .* HBx;
%             physClass.frhs{1}(:,:,3) = physClass.frhs{1}(:,:,3) + 2 * fphys{1}(:,:,6);
        end
        
    end
    
    methods( Access=protected )
        %> @brief Function to set up the matrix corresponding the Laplace operator and
        %> the first order operator related to nonhydrostatic pressure
        %> @details
        %> Function to set up the matrix corresponding the Laplace operator
        %> and the first order operator related to nonhydrostatic pressure
        %> @param[in] mesh The mesh object
        %> @param[in] element The studied element of a given mesh
        %> @param[in] edgeType the type of the edge of the studied element
        %> @param[out] PNPX the first order derivative of p with respect to x
        %> @param[out] SecondOrderTerm the discrete matrix of the Laplace operator regards to p
        [ PNPX, SecondOrderTerm ] = matAssembleCharacteristicMatrix( obj, mesh, element, edgeType)
        %> @brief Function to aspect the influence scope of a given study interpolation point
        %> @details
        %> Function to aspect the influence scope of a given study interpolation point,
        %> this could help to know advance the structure of the final stiff matrix, and
        %> helps to adopt openmp for accelaration
        %> @param[in] K The cell number of the studied mesh
        %> @param[in] Np The number of interpolation points per cell
        matAssemblePointToCellInformation( obj, K, Np )
        %> @brief Function to calculate the physical variable related partial derivative in a upwind manner
        %> @details
        %> Function to calculate the physical variable related partial derivative in a upwind manner
        %> @param[in] mesh The mesh object
        %> @param[in] fphys The fphys field
        %> @param[in] physClass The physical object set up
        %> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
        %> @param[out] UpwindedTermX The partial derivative with respect to x in upwind manner
        [ UpwindedTermX ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, termX)
        %> @brief Function to calculate the physical variable related partial derivative in a downwind manner
        %> @details
        %> Function to calculate the physical variable related partial derivative in a downwind manner
        %> @param[in] mesh The mesh object
        %> @param[in] fphys The fphys field
        %> @param[in] physClass The physical object set up
        %> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
        %> @param[out] DownwindedTermX The partial derivative with respect to x in downwind manner
        [ DownwindedTermX ] = matCalculateDownwindedFphysDerivative(obj, mesh, fphys, physClass, termX)
        %> @brief Function to calculate volume integral of the physical variable related partial derivative
        %> @details
        %> Function to calculate volume integral of the physical variable related partial derivative
        %> @param[in] mesh The mesh object
        %> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
        %> @param[out] VolumeX The volume integral of partial derivative with respect to x
        [ VolumeX ] = matVolumeIntegral( obj, mesh, termX )
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
        %> @brief Function to calculate the auxialary variable related partial derivative for IPDG
        %> @details
        %> Function to calculate the auxialary variable related partial derivative for IPDG
        %> @param[in] mesh The mesh object
        %> @param[in] BoundaryEdge The boundary edge object
        %> @param[in] InnerEdge The inner edge object
        %> @param[in] variable The primitive variable related to the auxialary variable
        %> @param[out] IPDGTermX The partial derivative with respect to x calculated with IPDG
        [ IPDGTermX ] = matCalculateIPDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, variable)
    end
    
end

