function GetCellSize( obj, mesh )
% Hcell = getStdCell( mesh.cell.N, enumStdCell.Line );
Hcell = StdLine(mesh.cell.N);
Vcell = StdLine(mesh.cell.Nz);
Vq = zeros( numel(Hcell.rq)*numel(Vcell.rq) ,( mesh.cell.N + 1 ) * ( mesh.cell.Nz + 1 ));
hq = repmat(Hcell.rq,numel(Vcell.rq),1);
vq = repmat(Vcell.rq,1,numel(Hcell.rq))';
hwq = repmat(Hcell.wq,numel(Vcell.wq),1);
vwq = repmat(Vcell.wq,1,numel(Hcell.wq))';
HInterp = ( Hcell.nodal_func(hq, zeros(size(hq)), zeros(size(hq))) );
VInterp = ( Vcell.nodal_func(vq(:), zeros(size(vq(:))), zeros(size(vq(:)))) );
for j = 1 : mesh.cell.Nz + 1
    for i = 1 : mesh.cell.N + 1
        Vq( :, (j - 1)*(mesh.cell.N + 1) + i ) = HInterp( :,i ) .* VInterp( :,j );
    end
end
wq = hwq.*vwq(:);
one = ones( obj.Nfp, obj.Ne );
obj.LAV = mxGetMeshIntegralValue( one, wq, obj.Js, Vq );
end