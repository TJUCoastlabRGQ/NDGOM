function [ Limiter ] = initSlopeLimiter2d(obj)

if obj.option.isKey('limiterType')
    type = obj.getOption('limiterType');
    [ Limiter ] = initSlopeLimiter( obj.mesh2d, type );
else % default non limiter
    Limiter = NdgNonLimiter( obj.mesh2d );
end

end

function [ limiter ] = initSlopeLimiter( mesh, type )
if ( type == enumLimiter.None )
    limiter = NdgNonLimiter( mesh );
elseif( type == enumLimiter.Vert )
    limiter = NdgVertLimiter2d( mesh );
elseif( type == enumLimiter.TVB )
    limiter = NdgTVB2d( mesh );
elseif( type == enumLimiter.BJ )
    limiter = NdgBJ2d( mesh );
end
end