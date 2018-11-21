classdef NdgVisAuxiStrongFormSolver < NdgQuadFreeStrongFormSolver & ...
        NdgVisAuxiAbstractSolver
    
    %NDGVISAUXISTRONGFORMSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgVisAuxiStrongFormSolver(phys)
            obj = obj@NdgQuadFreeStrongFormSolver(phys);
            obj = obj@NdgVisAuxiAbstractSolver(phys);
        end
        function evaluateViscosityRHS( obj, fphys )
            phys = obj.phys;
            phys.matUpdateVisParameter(fphys);
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion(m);
                [fm, fp] = phys.matEvaluateSurfaceValue(mesh, fphys{m}, phys.fext{m});
                varflux = obj.evaluateVarSurfNumFlux( fm, fp );
                deltaVarFlux =( -fm ) - ( -varflux );
                for i = 1:phys.Nvar
                    varFieldIndex = phys.varFieldIndex(i);
                    Epsilon = phys.epsilon{m}(:,:,varFieldIndex);
                    [qx, qy] = obj.evaluateAuxialaryVariable( m, fphys{m}(:,:,varFieldIndex), deltaVarFlux(:,:,i));
                    [Qx, Qy] = obj.evaluateVolumeFluxTerm(qx, qy, Epsilon);
                    [Numqx, Numqy] = obj.evaluateAuxialarySurfaceflux(mesh, Qx, Qy);
                    dqFlux = obj.nx{m} .* (Qx(mesh.eidM) - Numqx)  + obj.ny{m} .* (Qy(mesh.eidM) - Numqy);
                    [ phys.frhs{m}(:,:,i) ] = phys.frhs{m}(:,:,i)...
                        + obj.rx{m}.*( obj.Dr{m} * Qx ) ...
                        + obj.sx{m}.*( obj.Ds{m} * Qx ) ...
                        + obj.ry{m}.*( obj.Dr{m} * Qy ) ...
                        + obj.sy{m}.*( obj.Ds{m} * Qy ) ...
                        - ( obj.LIFT{m} * ( obj.Js{m} .* dqFlux ))./ obj.J{m};                    
                end
            end
        end

        %> @brief function used to calculate the auxialary variable, size(Np, K)
        %> @details
        %> @param[in] m Index of the studied mesh, size (TNfp,K,Nvar)
        %> @param[in] fphys The fphys value of the studied mesh of index varfieldIndex(i), size (Np,K)
        %> @retval[out] deltaVarflux the flux deviation of the surface flux term from numerical flux, size (TNfp,K)
        
        function [qx, qy] = evaluateAuxialaryVariable(obj,m, fphys, deltaVarflux)
            qx = obj.rx{m} .* (obj.Dr{m} * fphys)...
                + obj.sx{m} .* (obj.Ds{m} * fphys)...
                + obj.LIFT{m} * ( obj.nx{m} .* obj.Js{m} .* deltaVarflux )./obj.J{m};
            qy = obj.ry{m}.*(obj.Dr{m} * fphys)...
                + obj.sy{m}.* (obj.Ds{m} * fphys)...
                + obj.LIFT{m} * ( obj.ny{m} .* obj.Js{m} .* deltaVarflux )./obj.J{m};
        end
        
    end
    
end

