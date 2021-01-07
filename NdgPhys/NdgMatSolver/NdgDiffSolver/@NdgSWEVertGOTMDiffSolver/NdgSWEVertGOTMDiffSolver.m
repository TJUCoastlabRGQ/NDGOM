classdef NdgSWEVertGOTMDiffSolver < NdgVertDiffSolver
    %NDGVERTGOTMDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        h0b
    end
    
    properties
        GotmFile
    end
    
    methods
        function obj = NdgSWEVertGOTMDiffSolver( physClass )
            obj = obj@NdgVertDiffSolver( physClass );
            if  physClass.option.isKey('PhysicalBottomRoughnessLength')
                obj.h0b = physClass.getOption('PhysicalBottomRoughnessLength');
                fprintf('Value of the physical bottom roughness length is set to be: %f\n',obj.h0b);
            else
                obj.h0b = 0.05;
                fprintf('Value of the physical bottom roughness length is set to be the default value: %f\n',obj.h0b);
            end
            if  physClass.option.isKey('GOTMSetupFile')
                obj.GotmFile = physClass.getOption('GOTMSetupFile');
            else
                msg = 'Gotm setup file must be contained in the case folder, i.e. where the case begin';
                error(msg);
            end
            %Viscosity for GOTM turbulence model is initially set to be
            %zero
          obj.nv = zeros(size(physClass.meshUnion(1).x));
          obj.Prantl = physClass.Prantl;
          obj.matClearGlobalMemory;
        end
        
        function fphys = matUpdateImplicitVerticalDiffusion( obj, physClass, Height2d, Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage, Hu, Hv, time)
            obj.matUpdateViscosity( physClass, Height2d, Hu, Hv, dt, time);
            obj.matUpdatePenaltyParameter( physClass, obj.nv .* Height );
            fphys = obj.matCalculateImplicitRHS( physClass, obj.nv .* Height, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage);
        end
        
        function matClearGlobalMemory(obj)
            clear mxUpdateEddyViscosity;
        end
    end
    
    methods(Access = protected)
        function matUpdateViscosity(obj, physClass, H2d, Hu, Hv, dt, time)       
            [obj.nv, physClass.Cf{1}]  = mxUpdateEddyViscosity(physClass.mesh2d(1).cell.Np, physClass.mesh2d(1).K, physClass.meshUnion(1).cell.Np,...
                physClass.meshUnion(1).K, physClass.meshUnion(1).Nz, physClass.hcrit, physClass.finalTime, physClass.meshUnion(1).cell.VCV,...
                H2d, Hu, Hv, obj.GotmFile, dt, time, obj.h0b, physClass.SurfBoundNewmannDate(:,:,1), physClass.SurfBoundNewmannDate(:,:,2));        
        
        end
    end
    
end

