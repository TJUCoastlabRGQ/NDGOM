%> \brief determine the time interval
function [ dt ] = matUpdateTimeInterval( obj, fphys )
dt = nan;
for m = 1:obj.Nmesh
    dx = sqrt( obj.mesh2d(m).LAV );
    N = obj.mesh2d(m).cell.N;
    dtm = mxUpdateTimeInterval2d( obj.hcrit, ...
        obj.gra, N, dx, obj.mesh2d(m).status, fphys{m} );
    
    if ( dtm > 0 )
        dt = min(dt, dtm * obj.cfl);
    end
end

end