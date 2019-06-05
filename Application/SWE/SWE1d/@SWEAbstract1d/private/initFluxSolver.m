function [ numfluxSolver, surfluxSolver, volumefluxSolver ] = initFluxSolver( obj )
%INITNUMFLUXSOLVER Summary of this function goes here
%   Detailed explanation goes here

if obj.option.isKey('NumFluxType')
    if obj.option.isKey('nonhydrostaticType')
        if obj.getOption('NumFluxType') == enumSWENumFlux.HLL &&...
                obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
            % Numerical flux type pointed and non-hydrostatic type pointed
            numfluxSolver = SWEHLLNumFluxSolver1d( );
            surfluxSolver = SWEFaceFluxSolver1d( );
            % LF and Roe flux to be added
        elseif obj.getOption('NumFluxType') == enumSWENumFlux.HLL &&...
                obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Nonhydrostatic
            numfluxSolver = SWENonhydroHLLNumFluxSolver1d( );
            surfluxSolver = SWENonhydroFaceFluxSolver1d( );
            % LF and Roe flux to be added
        end    
    else  % No nonhydrostatic type pointed, and hydrostatic by default
        if obj.getOption('NumFluxType') == enumSWENumFlux.HLL 
            numfluxSolver = SWEHLLNumFluxSolver1d( );
            surfluxSolver = SWEFaceFluxSolver1d( );  
            % LF and Roe flux to be added
        end
    end
else  % No numerical flux type pointed
    if obj.option.isKey('nonhydrostaticType')
        if obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
            numfluxSolver = SWEHLLNumFluxSolver1d( );
            surfluxSolver = SWEFaceFluxSolver1d( );
        elseif obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Nonhydrostatic
            numfluxSolver = SWENonhydroHLLNumFluxSolver1d( );
            surfluxSolver = SWENonhydroFaceFluxSolver1d( );            
        end
    else % No nonhydrostatic type pointed, and hydrostatic by default
            numfluxSolver = SWEHLLNumFluxSolver1d( );
            surfluxSolver = SWEFaceFluxSolver1d( );          
    end
end
    parent = metaclass(obj);
    if strcmp(parent.SuperclassList.Name, 'SWEConventional1d')
        if obj.option.isKey('nonhydrostaticType')
            if obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
                volumefluxSolver = SWEVolumeFluxSolver1d();
            else
                volumefluxSolver = SWENonhydroVolumeFluxSolver1d();
            end
        else
            volumefluxSolver = SWEVolumeFluxSolver1d();
        end
    elseif strcmp(parent.SuperclassList.Name, 'SWEPreBlanaced1d') % SWEPreBlanaced1d
        if obj.option.isKey('nonhydrostaticType')
            if obj.getOption('nonhydrostaticType') == enumNonhydrostaticType.Hydrostatic
                volumefluxSolver = SWEPrebalanceVolumeFlux1d();
            else
                volumefluxSolver = SWENonhydroPrebalanceVolumeFluxSolver1d();
            end
        else
                volumefluxSolver = SWEPrebalanceVolumeFlux1d();      
        end       
    end
end

