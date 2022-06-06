classdef NdgVisCentralFluxSolver < NdgVisFluxAbstractSolver
    %NDGVISCENTRALFLUXSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        
        function obj = NdgVisCentralFluxSolver()
            obj = obj@NdgVisFluxAbstractSolver();
        end

        function  varflux = evaluateVarSurfNumFlux( obj, fm, fp )
            varflux = (fm + fp)/2;
        end
    
        function [Qx, Qy] = evaluateVolumeFluxTerm(obj, qx, qy, Epsilon)
                Qx = qx .* Epsilon;
                Qy = qy .* Epsilon;
        end
   
        function [Numqx, Numqy] = evaluateAuxialarySurfaceflux(obj, mesh, Qx, Qy)
            Numqx = ( Qx(mesh.eidM) + Qx(mesh.eidP) )/2;
            Numqy = ( Qy(mesh.eidM) + Qy(mesh.eidP) )/2;
        end
    end
    
end

