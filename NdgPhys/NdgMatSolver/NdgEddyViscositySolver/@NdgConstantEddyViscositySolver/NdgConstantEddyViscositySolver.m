classdef NdgConstantEddyViscositySolver < NdgAbstractEddyViscositySolver
    %	This solver stands for the constant eddy viscosity solver
    %   With this solver adopted, the constant eddy viscosity is set to be
    %   the user-difined constant value or the default 1.0e-4
    
    methods
        
        function obj = NdgConstantEddyViscositySolver( physClass )
                obj.matInitEddyViscosity(physClass,  physClass.mesh3d);
        end
        
        function  [EddyViscosity, Cf] = matUpdateEddyViscosity( obj, ~, ~, ~, ~, fphys, ~,~,~,~,Cf )
            EddyViscosity = fphys(:,:,5);
        end 
        
% matUpdateEddyViscosity( obj, ~, ~, ...
%         ~, ~, fphys{1}, ~ , ~, ~, ~, obj.Cf{1} );        
 
%      [fphys{1}(:,:,5), obj.Cf{1}] = obj.EddyViscositySolver.matUpdateEddyViscosity( obj, obj.mesh2d, ...
%         obj.meshUnion(1), fphys2d, fphys, dt , time, obj.WindTaux{1}, obj.WindTauy{1}, obj.Cf{1} );       
        
        
    end
    
    methods(Access = protected)
        function matInitEddyViscosity(obj, physClass, ~)
          if  physClass.option.isKey('ConstantEddyViscosityValue')
              value = physClass.getOption('ConstantEddyViscosityValue');
              fprintf('Value of the constant eddy viscosity is set to be: %f\n',value);
          else 
              value = 1.0e-4;
              fprintf('Value of the constant eddy viscosity is set to be the default value: %f\n',value);
          end
            for m = 1:physClass.Nmesh
                physClass.fphys{m}(:,:,5) = value;
            end
        end
    end
    
end