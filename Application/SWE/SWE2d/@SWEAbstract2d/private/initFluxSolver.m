function [ numfluxSolver, surfluxSolver, volumefluxSolver ] = initFluxSolver( obj )
%INITNUMFLUXSOLVER Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('NumFluxType')
    if obj.option.isKey('nonhydrostaticType')
        if obj.getOption('NumFluxType') == enumSWENumFlux.HLL &&...
                obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
            % Numerical flux type pointed and non-hydrostatic type pointed
            numfluxSolver = SWEHLLNumFluxSolver2d( );
            surfluxSolver = SWEFaceFluxSolver2d( );
            % LF and Roe flux to be added
        elseif obj.getOption('NumFluxType') == enumSWENumFlux.HLL &&...
                obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Nonhydrostatic
            numfluxSolver = SWENonhydroHLLNumFluxSolver2d( );
            surfluxSolver = SWENonhydroFaceFluxSolver2d( );
            % LF and Roe flux to be added
        end    
    else  % No nonhydrostatic type pointed, and hydrostatic by default
        if obj.getOption('NumFluxType') == enumSWENumFlux.HLL 
            numfluxSolver = SWEHLLNumFluxSolver2d( );
            surfluxSolver = SWEFaceFluxSolver2d( );  
            % LF and Roe flux to be added
        end
    end
else  % No numerical flux type pointed
    if obj.option.isKey('nonhydrostaticType')
        if obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
            numfluxSolver = SWEHLLNumFluxSolver2d( );
            surfluxSolver = SWEFaceFluxSolver2d( );
        elseif obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Nonhydrostatic
            numfluxSolver = SWENonhydroHLLNumFluxSolver2d( );
            surfluxSolver = SWENonhydroFaceFluxSolver2d( );            
        end
    else % No nonhydrostatic type pointed, and hydrostatic by default
            numfluxSolver = SWEHLLNumFluxSolver2d( );
            surfluxSolver = SWEFaceFluxSolver2d( );          
    end
end
    parent = metaclass(obj);
    if strcmp(parent.SuperclassList.Name, 'SWEConventional2d')
        if obj.option.isKey('nonhydrostaticType')
            if obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
                volumefluxSolver = SWEVolumeFluxSolver2d();
            else
                volumefluxSolver = SWENonhydroVolumeFluxSolver2d();
            end
        else
            volumefluxSolver = SWEVolumeFluxSolver2d();
        end
    elseif strcmp(parent.SuperclassList.Name, 'SWEPreBlanaced2d') % SWEPreBlanaced2d
        if obj.option.isKey('nonhydrostaticType')
            if obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
                volumefluxSolver = SWEPrebalanceVolumeFlux2d();
            else
                volumefluxSolver = SWENonhydroPrebalanceVolumeFluxSolver2d();
            end
        else
                volumefluxSolver = SWEPrebalanceVolumeFlux2d();      
        end       
    end
end

