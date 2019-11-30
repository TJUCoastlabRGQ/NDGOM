classdef VtkOutput2d < AbstractVtkOutput
    
methods (Access = public)
    function obj = VtkOutput2d( mesh, casename, Nfield, dt, varIndex )
        obj = obj@AbstractVtkOutput( mesh, casename, Nfield, dt, varIndex );
    end

    initFromMesh( obj, mesh );
    
    function closeOutputFile( obj )
    end
    
    % read output
    readOutputResult( obj, timeStep )
end
end