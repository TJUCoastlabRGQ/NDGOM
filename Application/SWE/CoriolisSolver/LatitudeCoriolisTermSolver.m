classdef LatitudeCoriolisTermSolver < AbstractCoriolisTermSolver
    
    properties (Constant)
        omega = 7.29e-5
    end
    
    properties (SetAccess = protected)
        LatField
        f
    end
    
    methods ( Access = public )
        function obj = LatitudeCoriolisTermSolver(phys, filename)
            obj = obj@AbstractCoriolisTermSolver( phys.meshUnion );
            LatVert = deg2rad( load(filename) );
            obj.LatField = cell( phys.Nmesh, 1 );
            obj.f = cell( phys.Nmesh, 1 );
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                obj.LatField{m} = mesh.proj_vert2node( LatVert );
                obj.f{m} = 2 * obj.omega * sin( obj.LatField{m} );
            end
        end
        
        function evaluateCoriolisTermRHS( obj, physClass, fphys )
            for m = 1:physClass.Nmesh
                mesh = physClass.meshUnion(m);
                ind = (mesh.status == int8(enumSWERegion.Wet));
                % frhs = frhs + f * hv
                physClass.frhs{m}(:, ind, obj.huIndex) = physClass.frhs{m}(:, ind, obj.huIndex)...
                    + obj.f{m}(:, ind) .* fphys{m}(:, ind, obj.hvIndex);
                
                % frhs = frhs - f * hu
                physClass.frhs{m}(:, ind, obj.hvIndex) = physClass.frhs{m}(:, ind, obj.hvIndex)...
                    - obj.f{m}(:, ind) .* fphys{m}(:, ind, obj.huIndex);
                
            end%for
        end
        
    end%methods
end
