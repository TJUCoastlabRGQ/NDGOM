classdef AbstractVtkOutput < AbstractOutputFile
    properties (SetAccess = protected)
        %> local node connection
        SEToV
        %> num of subcell
        Ncell
        %> num of points in each cell
        Np
        %>
        CellVertList
        %> num of connection
        Ncon
        %> points coordinate
        Points
        %> num of points
        Npoint
        %> cell type
        ctype
    end
    
    properties
        varName
        
        folderName
    end
    
    methods ( Access = public )
        function obj = AbstractVtkOutput( mesh, fieldName, casename, Nfield, dt, varIndex )
            obj = obj@AbstractOutputFile( mesh, casename, Nfield, dt, varIndex );
            obj.varName = cell(numel(varIndex),1);
            for i = 1:numel(varIndex)
                obj.varName{i} = fieldName{varIndex(i)};
            end
            if (mesh.type == enumMeshDim.One)
                obj.folderName = [casename,'/1d/'];
            elseif (mesh.type == enumMeshDim.Two)
                obj.folderName = [casename,'/2d/'];
            elseif (mesh.type == enumMeshDim.Three)
                obj.folderName = [casename,'/3d/'];
            end            
            
        end
        
        
    end
    
end