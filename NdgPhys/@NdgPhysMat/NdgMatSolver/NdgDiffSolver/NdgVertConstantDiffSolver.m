classdef NdgVertConstantDiffSolver < NdgVertDiffSolver
    %NDGVERTCONSTANTDIFFSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgVertConstantDiffSolver( physClass )
          if  physClass.option.isKey('ConstantVerticalEddyViscosityValue')
              value = physClass.getOption('ConstantVerticalEddyViscosityValue');
              fprintf('Value of the constant vertical eddy viscosity is set to be: %f\n',value);
          else 
              value = 1.0e-4;
              fprintf('Value of the constant vertical eddy viscosity is set to be the default value: %f\n',value);
          end
          obj.nv = value * ones(size(physClass.meshUnion(1).x));
          obj.Prantl = physClass.Prantl;
        end
    end
    
    methods(Access = protected)
        function matUpdateViscosity(obj, ~, ~, ~, ~, ~, ~)
            %doing nothing
        end
    end
    
end

