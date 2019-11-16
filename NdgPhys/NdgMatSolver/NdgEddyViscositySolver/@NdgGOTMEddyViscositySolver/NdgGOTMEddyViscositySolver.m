classdef NdgGOTMEddyViscositySolver < NdgAbstractEddyViscositySolver
    %NDGGOTMEDDYVISCOSITYSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        h0b
    end
    
    methods
        function obj = NdgGOTMEddyViscositySolver( physClass )
            if  physClass.option.isKey('PhysicalBottomRoughnessLength')
                value = physClass.getOption('PhysicalBottomRoughnessLength');
                disp('Value of the physical bottom roughness length is set to be: %f\n',value);
            else
                value = 0.05;
                disp('Value of the physical bottom roughness length is set to be the default value: %f\n',value);
            end
            
            obj.matInitEddyViscosity(physClass, physClass.mesh2d, physClass.mesh3d, physClass.hcrit);
        end
        
        function  EddyViscosity = matUpdateEddyViscosity( obj, fphys2d, fphys, dt, time  )  %fphys2d
            %             EddyViscosity = fphys(:,:,5);
            EddyViscosity;
        end
    end
    
    methods(Access = protected)
        function matInitEddyViscosity(obj, physClass, mesh2d, mesh3d, hcrit)
            %We note that, at present only one mesh can be considered when
            %initialize the GOTM model 
            if  physClass.option.isKey('GOTMSetupFile')
                file = physClass.getOption('GOTMSetupFile');
            else
                msg = 'Gotm setup file must be contained in the case folder, i.e. where the case begin';
                error(msg);
            end
           
            mxEddyViscosityByGOTMInit(mesh2d(1).K, mesh2d(1).cell.Np, mesh3d(1).K, mesh3d(1).cell.Np,...
                mesh3d(1).Nz, hcrit, obj.h0b, physClass.finalTime,file);
        end
    end
    
end

