classdef NdgSWEHorizSmagrinskyDiffSolver < NdgSWEHorizDiffSolver
    %NDGSWEHORIZSMAGRINSKYDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        C
    end
    
    methods
        function obj = NdgSWEHorizSmagrinskyDiffSolver( physClass )
            if physClass.SmagorinskyConstant == NULL
                obj.C = 100;
            else
                obj.C = physClass.SmagorinskyConstant;
            end
            obj.nv = zeros(size(physClass.meshUnion(1).x));
            obj.Prantl = physClass.Prantl;
        end
    end
    
    methods( Access = protected )
        function matUpdateViscosity( obj, physClass, hu, hv, h )
            Area = reshape(physClass.meshUnion(1).LAV, physClass.meshUnion(1).Nz, physClass.mesh2d.K)./physClass.meshUnion(1).Nz;
            obj.nv = 0.5 * obj.C .* Area .*...
                sqrt( ( physClass.meshUnion(1).rx .* ( physClass.meshUnion(1).cell.Dr * ( hu./h )) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds * ( hu./h )) ).^2 + ...
                0.5 * ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * ( hu./h )) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * ( hu./h )) + ...
                physClass.meshUnion(1).rx .* ( physClass.meshUnion(1).cell.Dr * ( hv./h )) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds * ( hv./h )) ).^2 +...
               ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * ( hv./h )) + ...
               physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * ( hv./h ))).^2 );
        end
    end
    
end

