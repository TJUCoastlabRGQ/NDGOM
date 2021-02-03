function [ fvert, fvmin, fvmax, cvar ] = matEvaluateVertAverage( obj, fphys, fieldId )

Nv = obj.Nv;
Nmesh = obj.Nmesh;
fvert = cell(Nmesh, 1);
fvmin = cell(Nmesh, 1);
fvmax = cell(Nmesh, 1);

% calculate the cell averages
cvar = cell( Nmesh, 1 );
for m = 1:Nmesh
    for i = 1:numel(fieldId)
        cvar{m}(:,:,i) = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,fieldId(i)) );
        %cvar{m} = obj.meshUnion(m).cell_mean( fphys{m}(:,:,fieldId) );
    end
end

for m = 1:Nmesh
    for i = 1:numel(fieldId)
        [ fvert{m}(:,:,i), fvmin{m}(:,:,i), fvmax{m}(:,:,i) ] = mxEvaluateVertAverage( cvar{m}(:,:,i), Nv, obj.Nvc, obj.VToM, obj.VToK, obj.VToW );
    end
end

end% func