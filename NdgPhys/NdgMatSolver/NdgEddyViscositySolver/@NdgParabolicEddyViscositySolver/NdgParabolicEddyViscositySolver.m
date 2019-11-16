classdef NdgParabolicEddyViscositySolver < NdgAbstractEddyViscositySolver
    %NDGPARABOLICEDDYVISCOSITYSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        
        function obj = NdgParabolicEddyViscositySolver( physClass )
            err = 'This is not implemented yet, choose eddy viscosity type Constant or GOTM';
%                 obj.matInitEddyViscosity(physClass, physClass.mesh2d, physClass.mesh3d, physClass.hcrit);
        end
        
        function  EddyViscosity = matUpdateEddyViscosity( obj, fphys, ~, ~  )
%             EddyViscosity = fphys(:,:,5);
        end
    end
    
    methods(Access = protected)
        function matInitEddyViscosity(obj, physClass, ~, ~, ~)
%           if  obj.option.isKey('ConstantEddyViscosityValue')
%               value = obj.getOption('ConstantEddyViscosityValue');
%               disp('Value of the constant eddy viscosity is set to be: %f\n',value);
%           else 
%               value = 1.0e-4;
%               disp('Value of the constant eddy viscosity is set to be the default value: %f\n',value);
%           end
%             for m = 1:physClass.meshUnion.Nmesh
%                 physClass.fphys{m}(:,:,5) = value;
%             end
        end
    end
    
end

