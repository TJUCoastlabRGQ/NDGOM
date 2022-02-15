classdef NdgSWEVertGOTMDiffSolver < NdgVertDiffSolver
    %NDGVERTGOTMDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        z0b
        
        z0s
    end
    
    properties
        Tke
        
        Eps
        
        rhoIndex
        
        maxhu
    end
    
    properties
        GotmFile
    end
    
    methods
        function obj = NdgSWEVertGOTMDiffSolver( physClass )
            obj = obj@NdgVertDiffSolver( physClass );
            if  physClass.option.isKey('PhysicalBottomRoughnessLength')
                obj.z0b = physClass.getOption('PhysicalBottomRoughnessLength');
                fprintf('Value of the physical bottom roughness length is set to be: %f\n',obj.z0b);
            else
                obj.z0b = 0.0015;
                fprintf('Value of the physical bottom roughness length is set to be the default value: %f\n',obj.z0b);
            end
            if  physClass.option.isKey('PhysicalSurfaceRoughnessLength')
                obj.z0s = physClass.getOption('PhysicalSurfaceRoughnessLength');
                fprintf('Value of the physical surface roughness length is set to be: %f\n',obj.z0s);
            else
                obj.z0s = 0.02;
                fprintf('Value of the physical surface roughness length is set to be the default value: %f\n',obj.z0s);
            end
            if  physClass.option.isKey('GOTMSetupFile')
                obj.GotmFile = physClass.getOption('GOTMSetupFile');
            else
                msg = 'Gotm setup file must be contained in the case folder, i.e. where the case begin';
                error(msg);
            end
            
            for i = 1:physClass.Nfield
                if (strcmp(physClass.fieldName3d{i},'rho'))
                    obj.rhoIndex = i;
                end
            end
            %Viscosity for GOTM turbulence model is initially set to be
            %zero
            obj.Prantl = physClass.Prantl;
            obj.matClearGlobalMemory;
            obj.maxhu = [];
        end
        
        function Outfphys = matUpdateImplicitVerticalDiffusion( obj, physClass, Height2d, Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage, Hu, Hv, time, fphys)
            Outfphys = obj.matCalculateImplicitRHS( physClass, obj.nv ./ Height./Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage, fphys{1}(:,:,1:2), Height2d);
            obj.matUpdateViscosity( physClass, Height2d, Hu, Hv, ImplicitParameter * dt, fphys{1}(:,:,obj.rhoIndex));
        end
        
        
        
        function matClearGlobalMemory(obj)
            matClearGlobalMemory@NdgVertDiffSolver( obj );
            clear mxUpdateEddyViscosity;
        end
    end
    
    methods(Access = protected)
        
        function matUpdateViscosity(obj, physClass, H2d, Hu, Hv, dt, rho )
            [ obj.nv, physClass.Cf{1}, obj.Tke, obj.Eps ]  = mxUpdateEddyViscosity(physClass.mesh2d(1).cell.Np, physClass.mesh2d(1).K, physClass.meshUnion(1).cell.Np,...
                physClass.meshUnion(1).K, physClass.meshUnion(1).Nz, physClass.hcrit, physClass.meshUnion(1).cell.VCV,...
                H2d, Hu, Hv, obj.GotmFile, dt, physClass.SurfBoundNewmannDate(:,:,1), physClass.SurfBoundNewmannDate(:,:,2), rho, obj.z0s, obj.z0b, physClass.gra, physClass.rho0, physClass.meshUnion.mesh2d.J,...
                physClass.meshUnion.mesh2d.cell.wq, physClass.meshUnion.mesh2d.cell.Vq, physClass.meshUnion.mesh2d.LAV );
        end           
        
        function matUpdataNewmannBoundaryCondition( obj, physClass, fphys)
            VCV = physClass.meshUnion(1).cell.VCV;
            Nz = physClass.meshUnion(1).Nz;
            Hu = VCV * fphys{1}(:,Nz:Nz:end,1);
            Hv = VCV * fphys{1}(:,Nz:Nz:end,2);
            H  = VCV * fphys{1}(:,Nz:Nz:end,4);
            physClass.BotBoundNewmannDate(:,:,1) = physClass.Cf{1} .* sqrt( (Hu./H).^2 + ...
                (Hv./H).^2 ) .* ( Hu./H ) * (-1);
            physClass.BotBoundNewmannDate(:,:,2) = physClass.Cf{1} .* sqrt( (Hu./H).^2 + ...
                (Hv./H).^2 ) .* ( Hv./H ) * (-1);
        end
    end
    
end

