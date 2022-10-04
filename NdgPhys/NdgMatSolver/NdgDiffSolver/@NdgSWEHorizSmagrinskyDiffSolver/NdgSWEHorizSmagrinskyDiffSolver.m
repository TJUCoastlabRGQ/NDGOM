classdef NdgSWEHorizSmagrinskyDiffSolver < NdgSWEHorizDiffSolver
    %NDGSWEHORIZSMAGRINSKYDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        C
    end
    
    methods
        function obj = NdgSWEHorizSmagrinskyDiffSolver( physClass )
            obj = obj@NdgSWEHorizDiffSolver( physClass );  
            if isempty(physClass.SmagorinskyConstant)
                obj.C = 0.01;
            else
                obj.C = physClass.SmagorinskyConstant;
            end
            obj.nv = zeros(size(physClass.meshUnion(1).x));
            obj.Prantl = physClass.Prantl;
        end
        
        function matEvaluateDiffRHS(obj, physClass, fphys)
            obj.matUpdateViscosity(physClass, fphys{1}(:,:,1), fphys{1}(:,:,2), fphys{1}(:,:,4));
            matEvaluateDiffRHS@NdgSWEHorizDiffSolver( obj, physClass, fphys);
        end
    end
    
    methods( Access = protected )
        function matUpdateViscosity( obj, physClass, hu, hv, h )
            obj.nv = obj.C .* bsxfun(@times, physClass.meshUnion(1).LAV./physClass.meshUnion(1).Nz ,...
                sqrt( ( physClass.meshUnion(1).rx .* ( physClass.meshUnion(1).cell.Dr * ( hu./h )) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds * ( hu./h )) ).^2 + ...
                0.5 * ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * ( hu./h )) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * ( hu./h )) + ...
                physClass.meshUnion(1).rx .* ( physClass.meshUnion(1).cell.Dr * ( hv./h )) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds * ( hv./h )) ).^2 +...
               ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * ( hv./h )) + ...
               physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * ( hv./h ))).^2 ));
        end
    end
    
end

