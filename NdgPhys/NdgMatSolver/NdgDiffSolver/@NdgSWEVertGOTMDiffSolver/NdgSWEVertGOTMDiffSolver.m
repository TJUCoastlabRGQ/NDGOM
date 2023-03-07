classdef NdgSWEVertGOTMDiffSolver < NdgSWEVertDiffSolver
    %NDGVERTGOTMDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        z0b
        
        z0s
    end
    
    properties
        Tke
        
        Eps
        
        nvh
    end
    
    properties
        GotmFile
    end
    
    properties
        uo
        
        vo
    end
    
    properties
        T0 = 0
        
        S0 = 0
         
        alphaT = 0
        
        betaS = 0
        
        EosType
    end
    
    methods
        function obj = NdgSWEVertGOTMDiffSolver( physClass )
            obj = obj@NdgSWEVertDiffSolver( physClass );
            
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
            if physClass.option.isKey('EosType')
                obj.EosType = physClass.getOption('EosType');
                if obj.EosType == enumEOSType.Linear
                    obj.T0 = physClass.T0;
                    obj.S0 = physClass.S0;
                    obj.alphaT = physClass.alphaT;
                    obj.betaS = physClass.betaS;
                end
            else
                msg = 'Type of the equation of state must be set first';
                error(msg);                
            end
            %Viscosity for GOTM turbulence model is initially set to be
            %zero
            obj.matClearGlobalMemory;
            obj.uo = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
            obj.vo = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
            obj.nvh = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
            obj.nv = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
            obj.BotBoundaryTreatType = 'Implicit';
        end
        
        function Outfphys = matUpdateImplicitVerticalDiffusion( obj, physClass, Height2d, Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage, Hu, Hv, time, fphys)
%             obj.matUpdataNewmannBoundaryCondition( physClass, fphys );
            Outfphys = obj.matCalculateImplicitRHS( physClass, (obj.nv + 1.3e-6) ./ Height./Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage, fphys{1}(:,:,1:2), Height2d, Height);
            obj.matUpdateViscosity( physClass, Height2d, obj.uo, obj.vo, Outfphys(:,:,1), Outfphys(:,:,2), zeros(size(Outfphys(:,:,2))), zeros(size(Outfphys(:,:,2))), ImplicitParameter * dt, Height, fphys{1}(:,:,14), fphys{1}(:,:,15) );
            Outfphys(:,:,3) = physClass.meshUnion.GetVertexWeightData( obj.Tke .* Height );
%             obj.Tke .* Height;
            Outfphys(:,:,4) = physClass.meshUnion.GetVertexWeightData( obj.Eps .* Height );
%             Outfphys = obj.matFilterData( physClass, Outfphys );
            obj.uo = Outfphys(:,:,1)./Height;
            obj.vo = Outfphys(:,:,2)./Height;
        end
        
        function matClearGlobalMemory(obj)
            matClearGlobalMemory@NdgSWEVertDiffSolver( obj );
            clear mxUpdateEddyViscosity;
        end
    end
    
    methods(Access = protected)
        
        function matUpdateViscosity(obj, physClass, H2d, uo, vo, HuNew, HvNew, HT, HS, dt, h, Hk, Heps )
            
            [ obj.nv, physClass.Cf{1}, obj.Tke, obj.Eps, obj.nvh ]  = mxUpdateEddyViscosity(physClass.mesh2d(1).cell.Np, physClass.mesh2d(1).K, physClass.meshUnion(1).cell.Np,...
                physClass.meshUnion(1).K, physClass.meshUnion(1).Nz, physClass.hcrit, physClass.meshUnion(1).cell.VCV,...
                H2d, uo, vo, obj.GotmFile, dt, physClass.SurfBoundNewmannDate(:,:,1), physClass.SurfBoundNewmannDate(:,:,2), ...
                obj.z0s, obj.z0b, physClass.gra, physClass.rho0, physClass.meshUnion.mesh2d.J, physClass.meshUnion.mesh2d.cell.wq,...
                physClass.meshUnion.mesh2d.cell.Vq, physClass.meshUnion.mesh2d.LAV, HuNew./h, HvNew./h,  physClass.meshUnion.J, ...
                physClass.meshUnion.cell.wq, physClass.meshUnion.cell.Vq, physClass.meshUnion.LAV, HT, HS, obj.T0, obj.S0, ...
                obj.alphaT, obj.betaS, char(obj.EosType), Hk./h, Heps./h);
                        
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
        
        function fphys = matFilterData( obj, physClass, fphys )
            fphys = mxFilterData( physClass.mesh2d(1).cell.Np, physClass.mesh2d(1).K, ...
                physClass.meshUnion(1).Nz, physClass.meshUnion.cell.Nz, fphys );
        end
        
    end
    
end

