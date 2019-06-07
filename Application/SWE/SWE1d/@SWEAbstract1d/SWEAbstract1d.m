classdef SWEAbstract1d < NdgPhysMat
    
    properties( Abstract, Constant )
        %> wet/dry depth threshold
        hmin
        %> gravity acceleration
        gra
    end
    
    properties
        %> number of physical field: [ h, hu, z, eta ]
        Nfield = 6
        %> number of variable field
        Nvar = 2
        %> index of variable in physical field
        varFieldIndex = [ 1, 2 ]
        %> order of the field to be written in the output file     
        outputFieldOrder = [1 2]
    end
    
    properties
        %> gradient of bottom elevation
        zGrad
        %> solver for evaluating friction source term
        frictionSolver
        %> solver for evaluating numerical flux
        numfluxSolver
        % solver for local face flux
        surfluxSolver
        %> solver for volume flux
        volumefluxSolver        
        %> slope limiter
        limiterSolver
    end
    
    properties( Hidden )
        draw_handle
    end
    
    methods
        %> draw physical field
        draw( obj, varargin );
        
        %> initilize from options
        initPhysFromOptions( obj, mesh );
        
        %> impose boundary condition and evaluate cell boundary values
        [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, fM, fP, fext );
    end
    
    methods( Abstract, Hidden )
        %> evaluate volume flux term
        [ E ] = matEvaluateFlux( obj, mesh, fphys )
    end
    
    methods( Abstract, Access = protected ) % Abstract function, protected function
        
        %> evaluate topography source term
        matEvaluateTopographySourceTerm( obj, fphys )
        
        %> evaluate wet/dry states
        matUpdateWetDryState(obj, fphys)
        
        %> apply post-process
        [ fphys ] = matEvaluatePostFunc(obj, fphys);
    end
        
    methods( Access = protected )
        
        %> apply slope limiter on conservative variables
        [ fphys ] = matEvaluateLimiter( obj, fphys )
        
        %> determine time interval
        [ dt ] = matUpdateTimeInterval( obj, fphys )
        
        %> evaluate source term
        [ ] = matEvaluateSourceTerm( obj, fphys )
        
        outputObj = matInitOutput( obj, mesh ) 
    end
    
    methods( Hidden, Sealed )
        %> evaluate boundary values
        function [ fM, fP ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
            [ fM, fP ] = mxEvaluateSurfaceValue1d( obj.hmin, obj.gra, ...
                mesh.eidM, mesh.eidP, mesh.nx, mesh.eidtype, fphys, fext );
        end
        
        %> evaluate local boundary flux
        function [ fluxM ] = matEvaluateSurfFlux( obj, mesh, nx, fm )
%             [ fluxM ] = mxEvaluateSurfFlux1d( obj.hmin, obj.gra, nx, fm);
            [ fluxM ] = obj.surfluxSolver.evaluate( obj.hmin, obj.gra, nx, fm);
        end
        
        %> evaluate boundary numerical flux
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, fm, fp )
            [ fluxS ] = obj.numfluxSolver.evaluate( obj.hmin, obj.gra, nx, fm, fp );
        end
    end
    
    %     methods( Access = protected )
    %         %>
    %         function [ fphys ] = matEvaluatePostFunc(obj, fphys)
    %             obj.matUpdateWetDryState( fphys );
    %             for m = 1:obj.Nmesh
    %                 hc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,1) );
    %                 qxc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,2) );
    %                 fphys{m}(:,:,1:2) = mxEvaluatePostFunc1d( obj.hmin, fphys{m}, hc, qxc );
    %             end
    %             obj.matUpdateWetDryState( fphys );
    %         end
    %     end
end

