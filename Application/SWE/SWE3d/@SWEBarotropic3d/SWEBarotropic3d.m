classdef SWEBarotropic3d < SWEAbstract3d
    %SWEBAROTROPIC3D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    %     properties ( Constant )
    %         %> number of physical field
    %         Nfield2d = 7 %[eta HU HV H Z Zx Zy]
    %         %> num of 3d physical field
    %         Nfield3d = 11 % [ Hu, Hv, w, tau_x, tau_y, H, eta, Z, Zx, Zy, MUy ]
    %         %> number of variable field
    %         Nvar3d = 2
    %         %> index of variable in physical field
    %         varFieldIndex3d = [ 1, 2 ]
    %         %> number of variable field
    %         Nvar2d = 1
    %         %> index of variable in physical field
    %         varFieldIndex2d = 1
    %     end
    properties
        %> number of physical field
        Nfield2d = 7 %[eta HU HV H Z Zx Zy]
        %> num of 3d physical field
        Nfield = 11 % [ Hu, Hv, w, tau_x, tau_y, H, eta, Z, Zx, Zy, MUy ]
        %> number of variable field
        Nvar = 2
        %> index of variable in physical field
        varFieldIndex = [ 1, 2 ]
        %> number of variable field
        Nvar2d = 1
        %> index of variable in physical field
        varFieldIndex2d = 1
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
            [ fluxS ] = obj.numfluxSolver.evaluate( obj, obj.hmin, obj.gra, nx, ny, nz, fm, fp );
        end% func
    end
    
    methods ( Access = protected )
        matEvaluateRK45( obj );
    end
    
end

