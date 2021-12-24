classdef NdgSWEHorizConstantDiffSolver < NdgSWEHorizDiffSolver
    %NDGSWEHORIZDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgSWEHorizConstantDiffSolver( physClass )
          obj = obj@NdgSWEHorizDiffSolver( physClass );  
          if  physClass.option.isKey('ConstantHorizontalEddyViscosityValue')
              value = physClass.getOption('ConstantHorizontalEddyViscosityValue');
              fprintf('Value of the constant horizontal eddy viscosity is set to be: %f\n',value);
          else 
              value = 0.01;
              fprintf('Value of the constant horizontal eddy viscosity is set to be the default value: %f\n',value);
          end
            obj.nv = value * ones(size(physClass.meshUnion(1).x));
            obj.Prantl = physClass.Prantl;
        end
        
        function matUpdateViscosity( ~ ,~ , ~, ~, ~ )
            %doing nothing
        end        
        
    end
    
    methods(Access = protected)

    end
end    
