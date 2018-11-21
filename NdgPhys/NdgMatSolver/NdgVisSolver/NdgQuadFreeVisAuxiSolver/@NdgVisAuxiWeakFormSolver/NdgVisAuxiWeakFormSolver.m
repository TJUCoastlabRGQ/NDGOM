classdef NdgVisAuxiWeakFormSolver < NdgQuadFreeWeakFormSolver & ...
        NdgVisAuxiAbstractSolver
    %NDGVISAUXIWEAKFORMSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgVisAuxiWeakFormSolver(phys)
            obj = obj@NdgQuadFreeWeakFormSolver(phys);
            obj = obj@NdgVisAuxiAbstractSolver(phys);
        end
        
        function evaluateViscosityRHS( obj, fphys )
            phys = obj.phys;
            phys.matUpdateVisParameter(fphys);
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion(m);
                [fm, fp] = phys.matEvaluateSurfaceValue(mesh, fphys{m}, phys.fext{m});
                varflux = obj.evaluateVarSurfNumFlux( fm, fp );

                for i = 1:phys.Nvar
                    varFieldIndex = phys.varFieldIndex(i);
                    Epsilon = phys.epsilon{m}(:,:,varFieldIndex);
                    [qx, qy] = obj.evaluateAuxialaryVariable( m, fphys{m}(:,:,varFieldIndex), varflux(:,:,i));
                    [Qx, Qy] = obj.evaluateVolumeFluxTerm(qx, qy, Epsilon);
                    [Numqx, Numqy] = obj.evaluateAuxialarySurfaceflux(mesh, Qx, Qy);
                    dqFlux  =  obj.nx{m} .* ( - Numqx(:,:,i)) + obj.ny{m} .*( - Numqy(:,:,i));
                    [ phys.frhs{m}(:,:,i) ] = phys.frhs{m}(:,:,i)...
                        - obj.Dr{m}*( obj.rx{m} .* (obj.M{m}*Qx) ) ...
                        - obj.Ds{m}*( obj.sx{m} .* (obj.M{m}*Qx) ) ...
                        - obj.Dr{m}*( obj.ry{m} .* (obj.M{m}*Qy) ) ...
                        - obj.Ds{m}*( obj.sy{m} .* (obj.M{m}*Qy) ) ...
                        - ( obj.LIFT{m} * ( obj.Js{m} .* dqFlux ))./ obj.J{m};
                end
            end
        end
        
        %> @brief function used to calculate the auxialary variables, size(Np, K)
        %> @details
        %> @param[in] m Index of the studied mesh
        %> @param[in] fphys The fphys value of the studied mesh of index varfieldIndex(i), size (Np,K)
        %> @retval[out] varflux the flux term of the variable numbered varFieldIndex(i) perpendicular to the interface, size (TNfp,K)
        
        function [qx, qy] = evaluateAuxialaryVariable(obj,m, fphys, Varflux)
            
            qx = -obj.Dr{m} * (obj.rx{m}.* (obj.M{m}*fphys))...
                - obj.Ds{m} * (obj.sx{m} .* (obj.M{m}*fphys))...
                + obj.LIFT{m} * ( obj.nx{m} .* obj.Js{m} .* Varflux )./obj.J{m};
            qy = -obj.Dr{m} * (obj.ry{m} .* (obj.M{m}*fphys))...
                - obj.Ds{m} * (obj.sy{m} .* (obj.M{m}*fphys))...
                + obj.LIFT{m} * ( obj.ny{m} .* obj.Js{m} .* Varflux )./obj.J{m};
            
        end
    end
    
end

