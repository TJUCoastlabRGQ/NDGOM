classdef VtkOutput1d < AbstractVtkOutput
    %VTKOUTPUT1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
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

