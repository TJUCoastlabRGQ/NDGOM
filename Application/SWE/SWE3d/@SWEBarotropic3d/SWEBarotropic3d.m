classdef SWEBarotropic3d < SWEAbstract3d
    %SWEBAROTROPIC3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %> number of physical field
        Nfield2d = 6 %[H HU HV Z Zx Zy]
        %> num of 3d physical field
        Nfield = 12 % [ Hu, Hv, w, H, miu, Z, eta, Zx, Zy ]
        
        fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy','w', 'hw','hc'};
        %> number of variable field
        Nvar = 2
        %> index of variable in physical field
        varFieldIndex = [ 1, 2 ]
        %> number of variable field
        Nvar2d = 1
        %> index of variable in physical field
        varFieldIndex2d = 1
        
        fieldName2d = {'h','hU','hV'};
    end
    
    properties
        %> the 2d field to be put in the output file
        outputFieldOrder2d =  1
        %> the 3d field to be put in the output file
        outputFieldOrder3d = [1 2 3]
    end
    
    methods
        function fphys = matImposeLimiter(obj, fphys)
            fphys = obj.Limiter.matLimit( fphys, obj.varFieldIndex(1) );
            fphys = obj.Limiter.matLimit( fphys, obj.varFieldIndex(2) );
        end
    end
    
    methods( Hidden )
        [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys );
        
    end
    
    methods ( Hidden, Access = public ) % public function, not allow to inherit
        
        [ fluxM ] = matEvaluateSurfFlux( obj, mesh, nx, ny, nz, fm )
        
        %> evaluate boundary numerical flux
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fm, fp, edge )
            [ fluxS ] = obj.numfluxSolver.evaluate( obj.hcrit, obj.gra, nx, ny, fm, fp, mesh, edge );
        end% func
    end
    
    methods ( Access = protected )
        
        matEvaluateRK45( obj );
        
        function matEvaluatePostFunc(obj, fphys2d)
            obj.matUpdateWetDryState(fphys2d);
        end
        
    end
    
end

