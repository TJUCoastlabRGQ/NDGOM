classdef NdgGOTMEddyViscositySolver < NdgAbstractEddyViscositySolver
    %NDGGOTMEDDYVISCOSITYSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        h0b
    end
        
    properties
        GotmFile
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
            
            obj.matInitEddyViscosity(physClass, physClass.mesh3d);
        end
        
        function  [ EddyViscosity, DragCoefficient ] = matUpdateEddyViscosity( obj, physClass, mesh2d, mesh3d, fphys2d, fphys, dt, time, WindTaux, WindTauy)  
            
            [EddyViscosity, DragCoefficient]  = mxUpdateEddyViscosity(mesh2d(1).cell.Np, mesh2d(1).K, mesh3d(1).cell.Np,...
                mesh3d(1).K, mesh3d(1).Nz, physClass.hcrit, physClass.finalTime, mesh3d(1).cell.VCV,...
                fphys2d{1}(:,:,1), fphys{1}(:,:,1), fphys{1}(:,:,2), obj.GotmFile, dt, time, obj.h0b, WindTaux, WindTauy);
            
        end
    end
    
    methods(Access = protected)
        function matInitEddyViscosity(obj, physClass, mesh3d)
            %We note that, at present only one mesh can be considered when
            %initialize the GOTM model
            
            if  physClass.option.isKey('GOTMSetupFile')
                obj.GotmFile = physClass.getOption('GOTMSetupFile');
            else
                msg = 'Gotm setup file must be contained in the case folder, i.e. where the case begin';
                error(msg);
            end
            
        end
    end
    
end

