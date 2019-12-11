classdef SWEBarotropic3d < SWEAbstract3d
    %SWEBAROTROPIC3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %> number of physical field
        Nfield2d = 6 %[H HU HV Z Zx Zy]
        %> num of 3d physical field
        Nfield = 9 % [ Hu, Hv, w, H, miu, Z, eta, Zx, Zy ]
        
        fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy'};
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
        outputFieldOrder = [1 2 3]
    end
    
    methods( Hidden )
        [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys );
    end
    
    methods ( Hidden, Access = public ) % public function, not allow to inherit
        
        [ fluxM ] = matEvaluateSurfFlux( obj, mesh, nx, ny, nz, fm )
        
        %> evaluate boundary numerical flux
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, nz, fm, fp )
            [ fluxS ] = obj.numfluxSolver.evaluate( obj, obj.hcrit, obj.gra, nx, ny, nz, fm, fp );
        end% func
    end
    
    methods ( Access = protected )
        matEvaluateRK45( obj );
    end
    
end

