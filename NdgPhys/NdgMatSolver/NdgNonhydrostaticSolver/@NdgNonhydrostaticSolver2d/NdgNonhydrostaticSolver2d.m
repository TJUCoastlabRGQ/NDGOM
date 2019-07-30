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
        
        PNPY
        
        TempPNPY
        
        HBy
        
        H2By
        
        fhy
        
        hvy
    end
    
    methods
        function obj = NdgNonhydrostaticSolver2d(PhysClass)
            
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
        %> @param[in] BoundaryEdgeny Projection of the unit direction vector in y direction of the boundary edge
        %> @param[in] BoundaryEdgeFToE The face to element topological relation of the boundary edge
        %> @param[in] InnerEdgenx Projection of the unit direction vector in x direction of the inner edge
        %> @param[in] InnerEdgeny Projection of the unit direction vector in y direction of the inner edge
        %> @param[in] InnerEdgeFToE The face to element topological relation of the inner edge
        %> @param[in] BoundaryEdgeJs The facial Jacobian of the boundary edge
        %> @param[in] InnerEdgeJs The facial Jacobian of the inner edge
        %> @param[out] nx Projection of the unit direction vector in x direction of the studied edge
        %> @param[out] ny Projection of the unit direction vector in y direction of the studied edge
        %> @param[out] Js The facial Jacobian of the studied edge
        [nx, ny, Js] = matGetElementFaceNormalDirectionVector( obj, ele, adjacentEle,...
            face, BoundaryEdgeFToF, BoundaryEdgenx, BoundaryEdgeny,...
            BoundaryEdgeFToE, InnerEdgenx, InnerEdgeny, InnerEdgeFToE, BoundaryEdgeJs, InnerEdgeJs );
        %> @brief Function to get the upwinded numerical flux when calculating the physical variable 
        %> related first order derivative in upwind manner
        %> @details
        %> Function to get the upwinded numerical flux when calculating the physical variable related first
        %> order derivative in upwind manner, we point out that this flux is set to be zero at the wet-dry 
        %> interface, such that the complex partial- wet cell problem can be avoided
        %> @param[in] status The mesh cell status
        %> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
        %> @param[in] HUm The HU field of the local face
        %> @param[in] HVm The HV field of the local face
        %> @param[in] HUp The HU field of the adjacent face
        %> @param[in] HVp The HV field of the adjacent face   
        %> @param[in] vfmx local face value of the studied variable in x direction
        %> @param[in] vfpx adjacent face value of the studied variable in x direction
        %> @param[in] vfmy local face value of the studied variable in y direction
        %> @param[in] vfpy adjacent face value of the studied variable in y direction    
        %> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
        %> @param[in] Edgeny Projection of the unit direction vector in y direction of the studied edge object   
        %> @param[out] fluxX the upwind numerical flux in x direction
        %> @param[out] fluxY the upwind numerical flux in y direction        
        [ fluxX, fluxY ] = matEvaluateUpwindNumFlux( obj, status, EdgeFToE, ...
            HUm, HVm, HUp, HVp, vfmx, vfpx, vfmy, vfpy, Edgenx, Edgeny);
        %> @brief Function to get the downwinded numerical flux when calculating the physical variable 
        %> related first order derivative in downwind manner
        %> @details
        %> Function to get the downwinded numerical flux when calculating the physical variable related first
        %> order derivative in downwind manner, we point out that this flux is set to be zero at the wet-dry 
        %> interface, such that the complex partial- wet cell problem can be avoided
        %> @param[in] status The mesh cell status
        %> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
        %> @param[in] HUm The HU field of the local face
        %> @param[in] HVm The HV field of the local face
        %> @param[in] HUp The HU field of the adjacent face
        %> @param[in] HVp The HV field of the adjacent face   
        %> @param[in] vfmx local face value of the studied variable in x direction
        %> @param[in] vfpx adjacent face value of the studied variable in x direction
        %> @param[in] vfmy local face value of the studied variable in y direction
        %> @param[in] vfpy adjacent face value of the studied variable in y direction    
        %> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
        %> @param[in] Edgeny Projection of the unit direction vector in y direction of the studied edge object   
        %> @param[out] fluxX the downwind numerical flux in x direction
        %> @param[out] fluxY the downwind numerical flux in y direction        
        [ fluxX, fluxY ] = matEvaluateDownwindNumFlux( obj, status, EdgeFToE, ...
            HUm, HVm, HUp, HVp, vfmx, vfpx, vfmy, vfpy, Edgenx, Edgeny);        
        %> @brief Function to get the local flux when calculating the physical variable related first order derivative
        %> @details
        %> Function to get the local flux when calculating the physical variable related first order derivative,
        %> we point out that this flux is set to be zero at the wet-dry interface, such that the complex partial-wet
        %> cell problem can be avoided
        %> @param[in] status The mesh cell status
        %> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
        %> @param[in] vfx local or adjacent face value of the studied variable in x direction
        %> @param[in] vfy local or adjacent face value of the studied variable in y direction
        %> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
        %> @param[in] Edgeny Projection of the unit direction vector in y direction of the studied edge object   
        %> @param[out] fluxMPx the local or adjacent flux in x direction
        %> @param[out] fluxMPy the local or adjacent flux in y direction            
        [ fluxMPx, fluxMPy ] = matEvaluateSurfFlux( obj, status, EdgeFToE, vfmx, vfmy, Edgenx, Edgeny);
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
        %> @param[out] PNPY the first order derivative of p with respect to y
        %> @param[out] SecondOrderTerm the discrete matrix of the Laplace operator regards to p
        [ PNPX, PNPY, SecondOrderTerm ] = matAssembleCharacteristicMatrix( mesh, element, edgeType)
        %> @brief Function to calculate the physical variable related partial derivative in a upwind manner
        %> @details
        %> Function to calculate the physical variable related partial derivative in a upwind manner
        %> @param[in] mesh The mesh object
        %> @param[in] fphys The fphys field
        %> @param[in] physClass The physical object set up
        %> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
        %> @param[in] termY The varialbe used to calculate the partial derivative with respect to y
        %> @param[out] UpwindedTermX The partial derivative with respect to x in upwind manner
        %> @param[out] UpwindedTermY The partial derivative with respect to y in upwind manner
        [ UpwindedTermX, UpwindedTermY ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, termX, termY)
        %> @brief Function to calculate the physical variable related partial derivative in a downwind manner
        %> @details
        %> Function to calculate the physical variable related partial derivative in a downwind manner
        %> @param[in] mesh The mesh object
        %> @param[in] fphys The fphys field
        %> @param[in] physClass The physical object set up
        %> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
        %> @param[in] termY The varialbe used to calculate the partial derivative with respect to y
        %> @param[out] DownwindedTermX The partial derivative with respect to x in downwind manner
        %> @param[out] DownwindedTermY The partial derivative with respect to y in downwind manner
        [ DownwindedTermX, DownwindedTermY ] = matCalculateDownwindedFphysDerivative(obj, mesh, fphys, physClass, termX, termY)
        %> @brief Function to calculate volume integral of the physical variable related partial derivative
        %> @details
        %> Function to calculate volume integral of the physical variable related partial derivative
        %> @param[in] mesh The mesh object
        %> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
        %> @param[in] termY The varialbe used to calculate the partial derivative with respect to y
        %> @param[out] VolumeX The volume integral of partial derivative with respect to x
        %> @param[out] VolumeY The volume integral of partial derivative with respect to y
        [VolumeX, VolumeY] = matVolumeIntegral( obj, mesh, termX, termY)
        %> @brief Function to calculate the auxialary variable related partial derivative for IPDG
        %> @details
        %> Function to calculate the auxialary variable related partial derivative for IPDG
        %> @param[in] mesh The mesh object
        %> @param[in] BoundaryEdge The boundary edge object
        %> @param[in] InnerEdge The inner edge object
        %> @param[in] variable The primitive variable related to the auxialary variable
        %> @param[out] IPDGTermX The partial derivative with respect to x calculated with IPDG
        %> @param[out] IPDGTermY The partial derivative with respect to y calculated with IPDG
        [IPDGTermX, IPDGTermY] = matCalculateIPDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, variable)
    end
end

