%> \brief
function obj = assembleEdgeConnect( obj, meshUnion )

mesh = obj.mesh;
K = mesh.K;
Nface = mesh.cell.Nface;
% IND has 3 rows, the first row is vert indice, 2nd row is element 
% indices, 3rd row is face indices.
ind = zeros( 3, K*Nface ); 

sk = 1;
for k = 1:K
    for f = 1:Nface
        locVertId = mesh.cell.FToV(:, f);
        ind( 1, sk ) = mesh.EToV(locVertId, k);
        ind( 2, sk ) = k;
        ind( 3, sk ) = f;
        sk = sk + 1;
    end
end
ind = ind';
[ ind ] = sortrows( ind );
ind = ind';

Nedge = 0;
flag = false( 1, K*Nface );
for i = 1:K*Nface
    if i == 1
        if any( ind(1,i+1) ~= ind(1,i) )
            Nedge = Nedge + 1;
            flag( i ) = 1;
        end
    elseif i == K*Nface
        if any( ind(1,i-1) ~= ind(1,i) )
            Nedge = Nedge + 1;
            flag( i ) = 1;
        end
    else
        if any( ind(1,i+1) ~= ind(1,i) ) && any( ind(1,i-1) ~= ind(1,i) )
            Nedge = Nedge + 1;
            flag( i ) = 1;
        end
    end
end

FToE = zeros( 2, Nedge);
FToF = zeros( 2, Nedge);
FToV = zeros( max( mesh.cell.Nfv ), Nedge );
FToM = ones( 2, Nedge ) * obj.FToM;

sk = 1;
for i = 1:K*Nface
    if( flag(i) )
        k = ind(2, i);
        f = ind(3, i);
        FToE(1, sk) = k;
        FToF(1, sk) = f;
        FToE(2, sk) = mesh.EToE(f, k);
        FToF(2, sk) = mesh.EToF(f, k);
        FToV(1, sk) = ind(1, i);
        FToM(2, sk) = mesh.EToM(f, k);
        sk = sk + 1;
    end
end

obj.Ne = Nedge;
obj.FToE = FToE;
obj.FToF = FToF;
obj.FToV = FToV;
obj.FToM = FToM;

end% func

% function vert = getFaceVertIndex( mesh, k, f )
% locVertId = mesh.cell.FToV(:, f);
% vert = sort( mesh.EToV(locVertId, k) );
% end