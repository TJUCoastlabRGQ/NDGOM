classdef NdgVertLimiter3d < NdgVertLimiter
    %NDGVERTLIMITER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        BoundaryEdge
        BottomBoundaryEdge
        SurfaceBoundaryEdge
    end
    
    methods
        
        function obj = NdgVertLimiter3d( mesh )
            obj = obj@NdgVertLimiter( mesh );
            warning('off');
            obj.BoundaryEdge = struct(mesh.BoundaryEdge);
            obj.BottomBoundaryEdge = struct(mesh.BottomBoundaryEdge);
            obj.SurfaceBoundaryEdge = struct(mesh.SurfaceBoundaryEdge);
            warning('on');
        end
        
        function fphys = matLimit( obj, fphys, fieldId )
            % calculate the cell averages
            Nmesh = obj.Nmesh;
            cvar = cell( Nmesh, 1 );
            for m = 1:Nmesh
                cvar{m} = obj.meshUnion(m).GetMeshAverageValue(...
                    fphys{m}(:,:,fieldId) );
                fphys{m}(:,:,fieldId) = mxVertLimit3d( fphys{m}(:,:,fieldId), cvar{m},...
                    obj.meshUnion(m).EToE, inv(obj.meshUnion(m).cell.V1d),...
                    obj.meshUnion(m).cell.Npz, obj.meshUnion(m).cell.Nph, ...
                    obj.meshUnion(m).cell.V1d);
            end
        end
        
        function fphys = matLimitNew( obj, physClass, fphys )
            
              fphys{1}(:,:,physClass.varFieldIndex) = ...
                mxVertLimit3dNew(fphys{1}, physClass.varFieldIndex, obj.meshUnion.cell.Np, ...
                obj.meshUnion.K, obj.meshUnion.J, obj.meshUnion.cell.wq, obj.meshUnion.cell.Vq, physClass.fext3d{1},...
                physClass.gra, physClass.hcrit, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                int8(physClass.meshUnion.BottomBoundaryEdge.ftype), ...
                int8(physClass.meshUnion.SurfaceBoundaryEdge.ftype), obj.BoundaryEdge, ...
                obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, obj.Nv, obj.Nvc, obj.VToK,...
                obj.meshUnion.mesh2d.cell.Fmask, obj.meshUnion.mesh2d.cell.Np,...
                obj.meshUnion.cell.Nz, obj.meshUnion.mesh2d.cell.Nv, obj.meshUnion.LAV,...
                obj.meshUnion.EToV);
        end
    end
    
end

