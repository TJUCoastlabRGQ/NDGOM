%> \brief determine the time interval
function [ dt ] = matUpdateTimeInterval( obj, fphys )
dt = nan;
for m = 1:obj.Nmesh
    dx = sqrt( obj.mesh2d(m).LAV );
    N = obj.mesh2d(m).cell.N;
    %     dtm = mxUpdateTimeInterval3d( obj.hcrit, ...
    %         obj.gra, N, dx, obj.mesh2d(m).status, fphys{m} );
    
    dtm = mxUpdateNonhydroTimeInterval3d(obj.gra, N, obj.meshUnion.cell.Nz, obj.meshUnion.cell.Np, dx, obj.mesh2d(m).status, fphys{m}(:,:,1),...
        fphys{m}(:,:,2), fphys{m}(:,:,4), fphys{m}(:,:,11), obj.meshUnion.Nz, obj.mesh2d(m).K );
    
    if ( dtm > 0 )
        dt = min(dt, dtm * obj.cfl);
    end
end

end
