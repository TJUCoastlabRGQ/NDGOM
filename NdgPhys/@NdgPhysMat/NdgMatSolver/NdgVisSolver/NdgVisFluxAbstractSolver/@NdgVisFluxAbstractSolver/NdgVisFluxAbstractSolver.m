classdef NdgVisFluxAbstractSolver < handle
    %NDGVISFLUXABSTRACTSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgVisFluxAbstractSolver()
        end
    end
    
    methods(Abstract)
        %> @brief function used to caculate the variable flux to compute the auxiallary flux
        %> @details
        %> @param[in] fm The surface value of the studied cell, size (TNfp,K,Nvar)
        %> @param[in] fp The surface value of the adjacent cell, size (TNfp,K,Nvar)
        %> @retval[out] varflux the flux term of the variable perpendicular to the interface, size (TNfp,K,Nvar)
        varflux = evaluateVarSurfNumFlux( obj, fm, fp );
        %> @brief function used to caculate the variable flux to compute the auxiallary volume flux
        %> @details
        %> @param[in] qx The auxiallary variable in the x direction
        %> @param[in] qy The auxiallary variable in the y direction
        %> @param[in] Epsilon The viscosity parameter of mesh m for varflux numbered varFieldIndex
        %> @retval[out] Qx, Qy the volume flux term of the auxiallary variable in each direction, size (Np,K)        
        [Qx, Qy] = evaluateVolumeFluxTerm(obj, qx, qy, Epsilon);
        %> @brief function used to caculate the auxiallary variable flux
        %> @details
        %> @param[in] mesh The mesh object        
        %> @param[in] Qx The volume flux term in x direction, size (TNfp,K)   
        %> @param[in] Qy The volume flux term in y direction, size (TNfp,K)   
        %> @retval[out] Numqx, Numqy the surface flux term of the auxiallary variable in each direction, size (TNfp,K)          
        [Numqx, Numqy] = evaluateAuxialarySurfaceflux(obj, mesh, Qx, Qy);
    end
    
end

