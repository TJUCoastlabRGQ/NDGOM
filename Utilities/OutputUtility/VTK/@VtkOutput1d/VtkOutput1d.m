classdef VtkOutput1d < AbstractVtkOutput
    %VTKOUTPUT1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
methods (Access = public)
    function obj = VtkOutput1d( casename, Nfield, dt )
        obj = obj@AbstractVtkOutput( casename, Nfield, dt );
    end

    initFromMesh( obj, mesh );
    
    function closeOutputFile( obj )
    end
    
    % read output
    readOutputResult( obj, timeStep )
end
    
end

