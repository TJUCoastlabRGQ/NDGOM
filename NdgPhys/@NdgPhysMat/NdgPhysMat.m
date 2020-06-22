%> @brief Abstract class for physical solver calculated by Matlab functions
%> 
%> The subclass of the NdgPhysMat use the Matlab or Mex functions to solve
%> the 
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef NdgPhysMat < NdgPhys

    properties
        %> cell array field for RHS
        frhs
    end
    
    properties (Abstract)
        %> order of the field to be written in the output file     
        outputFieldOrder
        
    end
    
    properties( SetAccess = protected )
        %> cell array for external value fields
        fext
        %> output netcdf file objects
        outputFile
        %> limiter object
        limiter
        %> 
        advectionSolver
        %> 
        viscositySolver
        %> 
        NonhydrostaticSolver
    end
    
    properties( SetAccess = protected )
        ftime
    end
    
    methods( Access = public )
        
        %> Evaluate the solver with the Matlab function
        function matSolve( obj )
            obj.matEvaluateTemporalDiscrete();
        end

        %> Initialize the solvers from input options and mesh
        initPhysFromOptions( obj, mesh );
    end
    
    methods
        %> @brief function for calculating the flux term
        %> @details
        %> Function to calculate the flux term of the equation.
        %> @param[in] mesh The ith mesh object
        %> @param[in] fphys The physical field on the ith mesh grid
        %> @retval[out] the flux term of the variable on each axis
        matEvaluateFlux( obj, mesh, fphys )
        %> @brief function for calculating the surface point values
        %> @details
        %> Function to calculate the flux term of the equation.
        %> @param[in] mesh The ith mesh object
        %> @param[in] fphys The physical field on the ith mesh grid
        %> @param[in] fext The external field on the ith mesh grid
        %> @retval[out] fm the variable in local element
        %> @retval[out] fp the variable on adjacent element
        [ fm, fp ] = matEvaluateSurfaceValue( obj, mesh, fphys, fext )
        %> @brief Evaluate the numerical flux
        %> @details
        %> Function calculates the flux division of the normal flux
        %> and the numerical flux term on mesh, which is obtained from
        %> \f$ \mathbf{F} \cdot \mathbf{n} - F^*  \f$
        %> @param [in] mesh The ith mesh object
        %> @param [in] fphys The physical field on the mesh grid
        %>
        [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, nz, fm, fp )
        %> @brief Evaluate the flux term on surface nodes
        %> @details
        %> Function calculates the flux division of the normal flux
        %> and the numerical flux term on mesh, which is obtained from
        %> \f$ \mathbf{F} \cdot \mathbf{n} - F^*  \f$
        %> @param [in] mesh The ith mesh object
        %> @param [in] fphys The physical field on the mesh grid
        %>
        [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, nz, fm )
        
        [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, ny, fM, fP, fext )
    end
    
    methods( Access = protected )
        %> @brief A warper functions for evaluating the RHS term
        %> @details
        %> For the 2d problem, the function should call function
        %> matEvaluateRHS2d; For the 3d problem, the function
        %> matEvaluateRHS3d should be called.
        matEvaluateRHS( obj, fphys )
        
        %> @brief Temporal discrete function
        %> @details
        %> The temporal discrete function will solve the function with the
        %> specific temporal discrete methods, which is defined in option
        %> by the name 'temporalDiscreteType'.
        function matEvaluateTemporalDiscrete( obj )            
            switch obj.getOption('temporalDiscreteType')
                case enumTemporalDiscrete.Euler
                    % call the Euler temporal discrete function
                    obj.matEvaluateEuler();
                case enumTemporalDiscrete.RK45
                    % call the SSP-RK45 temporal discrete function
                    obj.matEvaluateRK45();
                case enumTemporalDiscrete.RK22
                    obj.matEvaluateRK22();
                case enumTemporalDiscrete.RK33
                    obj.matEvaluateRK33();
                case enumTemporalDiscrete.SSPRK22
                    obj.matEvaluateSSPRK22();  
                case enumTemporalDiscrete.IMEXRK343 %This is implemented for barotropic swe3d only
                    obj.matEvaluateIMEXRK343();
                case enumTemporalDiscrete.IMEXRK222 %This is implemented for barotropic swe3d only
                    obj.matEvaluateIMEXRK222();    
                case enumTemporalDiscrete.IMEXRK111 %This is implemented for barotropic swe3d only
                    obj.matEvaluateIMEXRK111(); 
                case enumTemporalDiscrete.EXRK33 %This is implemented for swe2d only
                    obj.matEvaluateEXRK33();                    
                otherwise
                    msgID = [ mfilename, ':UnknownTemproalDicsreteType'];
                    msgtext = ['The temporal discrete type ', ...
                        obj.getOption('temporalDiscreteType') ,' is invalid.'];
                    throw( MException(msgID, msgtext) );
            end
        end% func
        
        %> An instantiation of the temporal-discrete function with the SSP-RK45 method
        matEvaluateRK45( obj );
        matEvaluateRK33( obj );
        matEvaluateRK22( obj );
        matEvaluateEXRK33( obj );
        %> An instantiation of the temporal-discrete with the Euler method
        matEvaluateEuler( obj );
        matEvaluateHeun( obj );
        matUpdateExternalField( obj, time, fphys )
        dt = matUpdateTimeInterval( obj, fphys )
        matEvaluateSourceTerm( obj, fphys )
        fphys = matEvaluateLimiter( obj, fphys )
        fphys = matEvaluatePostFunc( obj, fphys )
        
        outputObj = matInitOutput( obj )
        
        %> @brief
        matUpdateOutputResult( obj, time, step, fphys )
        
        %> @brief Update the final result of the physical field
        function matUpdateFinalResult( obj, time, fphys )
            for m = 1:obj.Nmesh
                obj.outputFile(m).outputFinalResult( ...
                    time, fphys{m}(:,:,obj.outputFile(m).varIndex) );
                obj.outputFile(m).closeOutputFile( ); %This function name can be changed again
            end
        end% function
    end
end
