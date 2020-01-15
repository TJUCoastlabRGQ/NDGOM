%> \brief determine the time interval
function [ dtm ] = matUpdateTimeInterval( obj, fphys2d )
% dt = nan;
for m = 1:obj.Nmesh
    dx = sqrt( obj.meshUnion(m).mesh2d.LAV );
    N = obj.meshUnion(m).mesh2d.cell.N;
    dtm = mxUpdateTimeInterval3d( obj.hmin, ...
        obj.gra, N, dx, obj.meshUnion(m).mesh2d.status, fphys2d{m} );
    
%     if ( dtm > 0 )
%         dt = min(dt, dtm * obj.cfl);
%     end
end

end

