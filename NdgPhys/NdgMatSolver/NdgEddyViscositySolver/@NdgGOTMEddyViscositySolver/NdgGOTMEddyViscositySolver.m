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
        
        function  EddyViscosity = matUpdateEddyViscosity( obj, physClass, mesh2d, mesh3d, fphys2d, fphys, dt  )  %fphys2d
            H2D = [1 2 3 4 5 6 7 8 9];
           Hu = 0.1:0.1:5.4;
           Hv = 0.1:0.1:5.4;
           dt = 0.03;
%            EddyViscosity = mxUpdateEddyViscosity(9, 1, 18, ...
%                3, 3, physClass.hcrit, mesh3d(1).cell.VCV, H2D, Hu, Hv, dt);           
%            EddyViscosity = mxUpdateEddyViscosity(mesh2d(1).cell.Np, mesh2d(1).K, mesh3d(1).cell.Np, ...
%                mesh3d(1).K, mesh3d(1).Nz, physClass.hcrit, mesh3d(1).cell.VCV, H2D, Hu, Hv, dt);
            EddyViscosity  = mxUpdateEddyViscosity(mesh2d(1).cell.Np, mesh2d(1).K, mesh3d(1).cell.Np,...
                mesh3d(1).K, mesh3d(1).Nz, physClass.hcrit, physClass.finalTime, mesh3d(1).cell.VCV,...
                fphys2d{1}(:,:,1), fphys{1}(:,:,1), fphys{1}(:,:,2), obj.GotmFile, dt, time, obj.h0b);
        end
    end
    
    methods(Access = protected)
        function matInitEddyViscosity(obj, physClass, mesh3d)
            %We note that, at present only one mesh can be considered when
            %initialize the GOTM model
            obj.TKE = zeros(mesh3d(1).Nz+1, physClass.mesh2d(1).K*physClass.mesh2d(1).cell.Np);
            obj.EPS = zeros(mesh3d(1).Nz+1, physClass.mesh2d(1).K*physClass.mesh2d(1).cell.Np);
            obj.L = zeros(mesh3d(1).Nz+1, physClass.mesh2d(1).K*physClass.mesh2d(1).cell.Np);
            obj.NUM = zeros(mesh3d(1).Nz+1, physClass.mesh2d(1).K*physClass.mesh2d(1).cell.Np);
            obj.NUH = zeros(mesh3d(1).Nz+1, physClass.mesh2d(1).K*physClass.mesh2d(1).cell.Np);
            
            physClass.fphys{1}(:,:,5) = 1.0e-8*ones(mesh3d(1).Np,mesh3d(1).K);
            
            if  physClass.option.isKey('GOTMSetupFile')
                obj.GotmFile = physClass.getOption('GOTMSetupFile');
            else
                msg = 'Gotm setup file must be contained in the case folder, i.e. where the case begin';
                error(msg);
            end

%             mxEddyViscosityByGOTMInit(mesh3d(1).Nz, file);            
            
%             mxEddyViscosityByGOTMInit(mesh2d(1).K, mesh2d(1).cell.Np, mesh3d(1).K, mesh3d(1).cell.Np,...
%                 mesh3d(1).Nz, hcrit, obj.h0b, physClass.finalTime,file);
        end
    end
    
end

