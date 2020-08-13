classdef NdgAbstractFilter < handle
    %NDGABSTRACTFILTER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %> the mesh object
        mesh
        %> number of mesh
        Nmesh
        %> The vandmonde matrix with high order basis set to be zero
        VF
        %> The smooth requirement parameter
        SR
        %> The diagnal filter strength matrix
        FS
    end
    
    properties
        %> the filter strength
        alpha = 0.1
        %> the filter order
        r = 2
    end
    
    methods
        
        function obj = NdgAbstractFilter( mesh )
            obj.mesh = mesh;
            obj.Nmesh = numel(mesh);
        end
        
        function fphys = matFilterSolution( obj, fphys, fieldId )
            for m = 1:obj.Nmesh
                fphys{m}(:,:,fieldId) = mxFilterSolution( obj.VF, obj.SR, obj.FS, obj.mesh.cell.wq, ...
                    obj.mesh.J, obj.mesh.cell.Vq, fphys{m}(:,:,fieldId), inv(obj.mesh.cell.V), obj.mesh.cell.V );
            end
        end
    end
    
    methods(Hidden, Abstract)
        
        matAssembleFilterVandMatirx( obj, mesh );
        
        matSetSmoothRequireParameter( obj, mesh );
        
        matSetFilterStrengthMatrix( obj, mesh );
        
    end
    
end

