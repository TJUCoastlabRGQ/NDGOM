function [ outputObj ] = matInitOutput( obj )
%INITOUTPUTFILE Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('nonhydrostaticType')
    switch obj.getOption('nonhydrostaticType')
        case enumNonhydrostaticType.Hydrostatic
            % Doing nothing
        case enumNonhydrostaticType.Nonhydrostatic
            % if Nonhydrostatic Solver included , the output variable is [ h, hu, hv, hw ]
            obj.Nvar = 4;
            obj.varFieldIndex = [1 2 3 6];
            obj.outputFieldOrder = [1 2 3 6 7];
    end% switch
end
outputObj = matInitOutput@NdgPhysMat(obj);

end

