classdef VtkOutput3d < VtkOutput2d
    methods (Access = public)
        function obj = VtkOutput3d( mesh, casename, Nfield, dt, varIndex3d )
            obj = obj@VtkOutput2d( mesh, casename, Nfield, dt, varIndex3d );
        end

        initFromMesh( obj, mesh );
        % drawResult( obj )
        
        readOutputResult( obj, timeStep )
    end
end