classdef NdgNonhydrostaticSolver1d < NdgAbstractNonhydrostaticSolver
    %NDGNONHYDROSTATICSOLVER1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function obj = NdgNonhydrostaticSolver1d(PhysClass)
            
            obj = obj@NdgAbstractNonhydrostaticSolver(PhysClass);
            
            obj.matSetInitializeCharacteristicMatrix( PhysClass.meshUnion(1) );
            
        end
        
    end
    
    methods ( Access = protected )
        %> @brief Function to set up the matrix corresponding the Laplace operator and
        %> the first order operator related to nonhydrostatic pressure
        %> @details
        %> Function to set up the matrix corresponding the Laplace operator
        %> and the first order operator related to nonhydrostatic pressure
        %> @param[in] mesh The mesh object
        matSetInitializeCharacteristicMatrix(obj, mesh)
        %> @brief Function to make the nonhydrostatic correction
        %> @details
        %> Function to make the nonhydrostatic correction
        %> @param[in] physClass The physical object set up
        %> @param[in] fphys The fphys field
        fphys = matNdgConservativeNonhydrostaticUpdata( obj, physClass, fphys )
        %> @brief Function to reassemble the characteristic matrix
        %> @details
        %> Function to reassemble the characteristic matrix with wetting and drying considered
        matReconstructStiffmatrixRelatedMatrix(obj)
        %> @brief Function to calculate the physical variable related partial derivative
        %> @details
        %> Function to calculate the physical variable related partial derivative
        %> @param[in] mesh The mesh object
        %> @param[in] fphys The fphys field
        %> @param[in] physClass The physical object set up
        matCalculateFphysDerivative( obj, mesh, fphys, physClass )
        %> @brief Function to get the unit vector and facial Jacobian at
        %> the studied face
        %> @details
        %> Function to get the unit vector and facial Jacobian at the studied face
        %> @param[in] ele The studied local element
        %> @param[in] adjacentEle The studied adjacent element
        %> @param[in] face The studied local face
        %> @param[in] BoundaryEdgeFToF The face to face topological relation of the boundary edge
        %> @param[in] BoundaryEdgenx Projection of the unit direction vector in x direction of the boundary edge
        %> @param[in] BoundaryEdgeFToE The face to element topological relation of the boundary edge
        %> @param[in] InnerEdgenx Projection of the unit direction vector in x direction of the inner edge
        %> @param[in] InnerEdgeFToE The face to element topological relation of the inner edge
        %> @param[in] BoundaryEdgeJs The facial Jacobian of the boundary edge
        %> @param[in] InnerEdgeJs The facial Jacobian of the inner edge
        %> @param[out] nx Projection of the unit direction vector in x direction of the studied edge
        %> @param[out] Js The facial Jacobian of the studied edge
        [nx, Js] = matGetElementFaceNormalDirectionVector( obj, ele, adjacentEle,...
            face, BoundaryEdgeFToF, BoundaryEdgenx, BoundaryEdgeFToE, InnerEdgenx,...
            InnerEdgeFToE, BoundaryEdgeJs, InnerEdgeJs );
        %> @brief Function to get the upwinded numerical flux when calculating the physical variable
        %> related first order derivative in upwind manner
        %> @details
        %> Function to get the upwinded numerical flux when calculating the physical variable related first
        %> order derivative in upwind manner, we point out that this flux is set to be zero at the wet-dry
        %> interface, such that the complex partial- wet cell problem can be avoided
        %> @param[in] status The mesh cell status
        %> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
        %> @param[in] HUm The HU field of the local face
        %> @param[in] HUp The HU field of the adjacent face
        %> @param[in] vfmx local face value of the studied variable in x direction
        %> @param[in] vfpx adjacent face value of the studied variable in x direction
        %> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
        %> @param[out] fluxX the upwind numerical flux in x direction
        [ fluxX ] = matEvaluateUpwindNumFlux( obj, status, EdgeFToE, ...
            HUm, HUp, vfmx, vfpx, Edgenx);
        %> @brief Function to get the downwinded numerical flux when calculating the physical variable
        %> related first order derivative in downwind manner
        %> @details
        %> Function to get the downwinded numerical flux when calculating the physical variable related first
        %> order derivative in downwind manner, we point out that this flux is set to be zero at the wet-dry
        %> interface, such that the complex partial- wet cell problem can be avoided
        %> @param[in] status The mesh cell status
        %> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
        %> @param[in] HUm The HU field of the local face
        %> @param[in] HUp The HU field of the adjacent face
        %> @param[in] vfmx local face value of the studied variable in x direction
        %> @param[in] vfpx adjacent face value of the studied variable in x direction
        %> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
        %> @param[out] fluxX the downwind numerical flux in x direction
        [ fluxX ] = matEvaluateDownwindNumFlux( obj, status, EdgeFToE, ...
            HUm, HUp, vfmx, vfpx, Edgenx);
        %> @brief Function to get the local flux when calculating the physical variable related first order derivative
        %> @details
        %> Function to get the local flux when calculating the physical variable related first order derivative,
        %> we point out that this flux is set to be zero at the wet-dry interface, such that the complex partial-wet
        %> cell problem can be avoided
        %> @param[in] status The mesh cell status
        %> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
        %> @param[in] vfx local or adjacent face value of the studied variable in x direction
        %> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
        %> @param[out] fluxMPx the local or adjacent flux in x direction
        [ fluxMPx ] = matEvaluateSurfFlux( obj, status, EdgeFToE, vfmx, Edgenx);
        %> @brief Function to get the inner edge that cells located besides are wet and non-wet cell
        %> @details Function to get the inner edge that cells located besides are wet and non-wet cell
        %> @param[in] status Status of the mesh cell
        %> @param[in] InnerEdgeFToE The face to element topological relation of the inner edge
        %> @param[out] faceflag the face status flag, 1 for wet-dry interface, 0 for wet-wet interface        
        faceflag = matGetWetDryFace( obj, status, InnerEdgeFToE)        
    end
    
    
    methods(  Access = protected, Abstract )
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
        [ PNPX, SecondOrderTerm ] = matAssembleCharacteristicMatrix( mesh, element, edgeType)
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

