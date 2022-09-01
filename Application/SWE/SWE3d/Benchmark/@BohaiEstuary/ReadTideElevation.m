function ReadTideElevation( obj )
%READTIDEELEVATION
Tempdata = xlsread( obj.tidalFile );

data = zeros(size(Tempdata,1),2 + (size(Tempdata,2) - 2)*2);

data(:,1) = Tempdata(:,1);
for i = 2:size(Tempdata,2)-2
    data(:,1+(i-2)*2+1) = Tempdata(:,i);
    data(:,1+(i-2)*2+2) = Tempdata(:,i);
end
data(:,end) = Tempdata(:,end);

TotalNb = 0;
Nb = cell( obj.Nmesh, 1 );
obj.OBEdgeIndex = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    edge = obj.mesh2d( m ).BoundaryEdge;
    obj.OBEdgeIndex{m} = find( edge.ftype == enumBoundaryCondition.ClampedDepth );
    Nb{m} = edge.Nfp * numel( obj.OBEdgeIndex{m} );
    TotalNb = TotalNb + Nb{m};
end
% Ntime = numel( data ) / TotalNb;
Ntime = size(data,1);

obj.Tide = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    Nbedge = numel( obj.OBEdgeIndex{m} );
    obj.Tide{m} = zeros( edge.Nfp, Nbedge, Ntime );
end

% counter = 1;
for t = 1:Ntime
    for m = 1:obj.Nmesh
        mesh = obj.mesh2d( m );
        edge = obj.mesh2d( m ).BoundaryEdge;
        
        ind = edge.FToN1 + repmat(edge.FToE(1, :) - 1, edge.Nfp, 1) * mesh.cell.Np;
        bot = obj.fphys2d{m}(:, :, 4);
        Index = ind(:,obj.OBEdgeIndex{m});
%         obj.fext{m}(:, :, 4) = bot( ind );
        
        Nbedge = numel( obj.OBEdgeIndex{m} );
        obj.Tide{m}(:, :, t) = reshape( ...
            data(t,:) - bot(Index(:))', edge.Nfp, Nbedge );
%         obj.Tide{m}(:, :, 1436) = reshape( ...
%             - bot(Index(:))', edge.Nfp, Nbedge );        
%         obj.Tide{m}(:, :, t) = reshape( ...
%             - bot(Index(:))', edge.Nfp, Nbedge );  
%         counter = counter + Nb{m};
    end
end

end
