classdef NdgSWEVertConstantDiffSolver < NdgVertDiffSolver
    %NDGVERTDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        
        function obj = NdgSWEVertConstantDiffSolver( physClass )
            obj = obj@NdgVertDiffSolver( physClass );
          if  physClass.option.isKey('ConstantVerticalEddyViscosityValue')
              value = physClass.getOption('ConstantVerticalEddyViscosityValue');
              fprintf('Value of the constant vertical eddy viscosity is set to be: %f\n',value);
          else 
              value = 0;
              fprintf('Value of the constant horizontal eddy viscosity is set to be the default value: %f\n',value);
          end
            obj.nv = value * ones(size(physClass.meshUnion(1).x));
            obj.Prantl = physClass.Prantl;            
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
        %         fphys = matUpdateImplicitVerticalDiffusion( obj, physClass, Height2d, Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage, Hu, Hv, time)
        function fphys = matUpdateImplicitVerticalDiffusion( obj, physClass, ~, Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage, ~, ~, ~, ~)
%             obj.matUpdatePenaltyParameter( physClass, obj.nv ./ Height.^2 );
            fphys = obj.matCalculateImplicitRHS( physClass, obj.nv ./ Height.^2, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage);
        end
    end
    
end

