classdef NdgSWEVertDiffSolver < handle
    %NDGSWEVERTDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        BoundaryEdgeType = 'Dirichlet'
    end
    
    properties (Constant)
        omega = 7.29e-5
    end
    
    properties
        nv
        
        %> the Prantal number for heat and sault diffusion
        Prantl
    end
    
    properties
        % How to treat the bottom Neumann boundary, Implicit or Explicit
        BotBoundaryTreatType = 'Implicit'
    end
    
    properties
        f0
    end
    
    methods
        function obj = NdgSWEVertDiffSolver(physClass)
            if physClass.option.isKey('CoriolisType') % the option exist
                switch physClass.getOption('CoriolisType')
                    case enumSWECoriolis.None
                        obj.f0 = zeros(size(physClass.meshUnion.x));
                    case enumSWECoriolis.Beta
                        if physClass.option.isKey('f0 for beta coriolis solver') && physClass.option.isKey('beta for beta coriolis solver')
                            f = physClass.getOption('f0 for beta coriolis solver');
                            beta = physClass.getOption('beta for beta coriolis solver');
                            obj.f0 = f + beta*physClass.meshUnion.y;
                        else
                            msgID = [ mfilename, ':Parameter required not supplied'];
                            msgtext = ['Parameter f0 and beta should be given in the set up part for the beta approximation coriolis solver'];
                            throw( MException(msgID, msgtext) );
                        end
                    case enumSWECoriolis.Latitude
                        if physClass.option.isKey('Latitude file')
                            filename = physClass.getOption('Latitude file');
                            LatVert = deg2rad( load(filename) );
                            LatField = physClass.meshUnion.proj_vert2node( LatVert );
                            obj.f0 = 2 * obj.omega * sin( LatField );
                        else
                            msgID = [ mfilename, ':Parameter required not supplied'];
                            msgtext = ['Latitude file should be given in the set up part for the latitude coriolis solver'];
                            throw( MException(msgID, msgtext) );
                        end
                    otherwise
                        msgID = [ mfilename, ':Wrong type coriolis solver'];
                        msgtext = ['The coriolis solver type is invalid.'];
                        throw( MException(msgID, msgtext) );
                end
                %doing nothing
            else % the option does not exist
                obj.f0 = zeros(size(physClass.meshUnion.x));
            end
            if physClass.option.isKey('BottomBoundaryEdgeType')
                obj.BoundaryEdgeType = char(physClass.getOption('BottomBoundaryEdgeType'));
            end
            fprintf('The bottom boundary condition for momentum is: %s\n',obj.BoundaryEdgeType);
            
            obj.matClearGlobalMemory;
        end
        
        function matClearGlobalMemory( obj )
            clear mxUpdateImplicitRHS;
            clear mxSparseVersionUpdateImplicitRHS;
        end
        
        %> @brief Calculating the right hand side corresponding to the vertical diffusion term and
        %> return the physical field with vertical diffusion considered
        %> @detail this function is used to calculate the right hand side corresponding to the vertical
        %> diffusion term and return the updated physical field at each Runge-Kutta time stage
        %> @param[in] physClass The physical solver establised
        %> @param[in] Height The water depth
        %> @param[in] ImplicitParameter The implicit parameter at the corresponding IMEXRK stage
        %> @param[in] dt The time step
        %> @param[out] fphys The physical field with vertical diffusion
        %> considered
        %> Input parameter changed on 20211231 to consider the hu and hv
        %> field, since we need it when we treat the bottom boundary implicitly.
        function  fphys = matUpdateImplicitVerticalDiffusion( obj, physClass, SystemRHS, ImplicitParameter, dt, intRK, Stage, huv3d, h2d )
            fphys = obj.matCalculateImplicitRHS( physClass, obj.nv, SystemRHS, ImplicitParameter, dt, intRK, Stage, huv3d, h2d);
        end
        
    end
end

