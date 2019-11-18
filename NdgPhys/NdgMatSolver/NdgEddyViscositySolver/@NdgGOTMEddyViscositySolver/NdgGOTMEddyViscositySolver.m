classdef NdgGOTMEddyViscositySolver < NdgAbstractEddyViscositySolver
    %NDGGOTMEDDYVISCOSITYSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
        h0b
    end
    
    methods
        function obj = NdgGOTMEddyViscositySolver( physClass )
            if  physClass.option.isKey('PhysicalBottomRoughnessLength')
                obj.h0b = physClass.getOption('PhysicalBottomRoughnessLength');
                fprintf('Value of the physical bottom roughness length is set to be: %f\n',obj.h0b);
            else
                obj.h0b = 0.05;
                fprintf('Value of the physical bottom roughness length is set to be the default value: %f\n',obj.h0b);
            end
            
            obj.matInitEddyViscosity(physClass, physClass.mesh2d, physClass.mesh3d, physClass.hcrit);
        end
        
        function  EddyViscosity = matUpdateEddyViscosity( obj, mesh3d, fphys2d, fphys, dt, time  )  %fphys2d
            H2D = [1 2 3 4 5 6 7 8 9];
           Hu = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
           Hv = [0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2];
           dt = 0.03;
           time = 5;
           EddyViscosity = mxUpdateEddyViscosity(mesh3d(1).cell.VCV, H2D, Hu, Hv, dt, time);
%             EddyViscosity = mxUpdateEddyViscosity(mesh3d(1).cell.VTV, fphys2d{1}(:,:,1), fphys{1}(:,:,1), fphys{1}(:,:,2),...
%                 dt, time);
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

