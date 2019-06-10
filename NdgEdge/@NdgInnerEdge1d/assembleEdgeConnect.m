%>
function obj = assembleEdgeConnect( obj, meshUnion )

mesh = obj.mesh;
K = mesh.K;
Nface = mesh.cell.Nface;
% IND has 3 rows, first row is vert indice, 2nd row is element
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
ind = sortrows( ind );
ind = ind';

Nedge = 0;
flag = false( 1, K*Nface );

% check first 2 faces
if ind( 1, 2 ) == ind( 1, 1 )
    Nedge = Nedge + 1;
    flag( 1 ) = 1;
end

for i = 2 : ( K*Nface - 1 )
    if ind( 1, i+1 ) == ind( 1, i ) ...
            && ind( 1, i-1 ) ~= ind( 1, i ) 
        Nedge = Nedge + 1;
        flag( i ) = 1;
    end
end

FToE = zeros(2, Nedge);
FToF = zeros(2, Nedge);
FToV = zeros( max( mesh.cell.Nfv ), Nedge );

sk = 1;
for i = 1:K*Nface
    if( flag(i) )
        FToE(1, sk) = ind(2, i);
        FToF(1, sk) = ind(3, i);
        FToE(2, sk) = ind(2, i+1);
        FToF(2, sk) = ind(3, i+1);
        FToV(1, sk) = ind(1, i);
        sk = sk + 1;
    end
end

obj.Ne = Nedge;
obj.FToE = FToE;
obj.FToF = FToF;
obj.FToV = FToV;

end% func

% function vert = getFaceVertIndex( mesh, k, f )
% locVertId = mesh.cell.FToV(:, f);
% vert = sort( mesh.EToV(locVertId, k) );
% %ind = sum( vert .* (mesh.Nv .^ ((numel(vert)-1):-1:0)') );
% end
