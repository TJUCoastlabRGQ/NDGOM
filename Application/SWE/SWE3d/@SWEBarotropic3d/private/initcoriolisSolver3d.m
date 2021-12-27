function coriolisSolver3d = initcoriolisSolver3d(obj)
if obj.option.isKey('CoriolisType') % the option exist
    switch obj.getOption('CoriolisType')
        case enumSWECoriolis.None
            coriolisSolver3d = NonCoriolisTermSolver();
        case enumSWECoriolis.Beta
            if obj.option.isKey('f0 for beta coriolis solver') && obj.option.isKey('beta for beta coriolis solver')
                f0 = obj.getOption('f0 for beta coriolis solver');
                beta = obj.getOption('beta for beta coriolis solver');
            else
                msgID = [ mfilename, ':Parameter required not supplied'];
                msgtext = ['Parameter f0 and beta should be given in the set up part for the beta approximation coriolis solver'];
                throw( MException(msgID, msgtext) );
            end
            coriolisSolver3d = BetaApproCoriolisTermSolver( obj, f0, beta );
        case enumSWECoriolis.Latitude
            if obj.option.isKey('Latitude file')
                filename = obj.getOption('Latitude file');
            else
                msgID = [ mfilename, ':Parameter required not supplied'];
                msgtext = ['Latitude file should be given in the set up part for the latitude coriolis solver'];
                throw( MException(msgID, msgtext) );
            end
            coriolisSolver3d = LatitudeCoriolisTermSolver(obj, filename);
        otherwise
            msgID = [ mfilename, ':Wrong type coriolis solver'];
            msgtext = ['The coriolis solver type is invalid.'];
            throw( MException(msgID, msgtext) );
    end
    %doing nothing
else % the option does not exist
    coriolisSolver3d = NonCoriolisTermSolver();
end
end