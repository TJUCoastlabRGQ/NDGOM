function [ outputObj ] = matInitOutput( obj, mesh )
%INITOUTPUTFILE Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('nonhydrostaticType')
    switch obj.getOption('nonhydrostaticType')
        case enumNonhydrostaticType.Hydrostatic
            % Doing nothing
        case enumNonhydrostaticType.Nonhydrostatic
            % if Nonhydrostatic Solver included , the output variable is [ h, hu, hw ]
            obj.Nvar = 3;
            obj.varFieldIndex = [1 2 5];
            obj.outputFieldOrder = [1 2 5 6];
    end% switch
end

outputObj = matInitOutput@NdgPhysMat(obj, mesh);

end

