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
                %The following way can be only used for constant Jacobi
                Data = obj.meshUnion.cell.V\fphys{m}(:,:,fieldId);
                data = Data(1,:) .* obj.meshUnion.cell.V(1) ;
                difference = data - cvar{m};
               fphys{m}(:,:,fieldId) = mxVertLimit3d( fphys{m}(:,:,fieldId), cvar{m},...
                   obj.meshUnion(m).EToE, inv(obj.meshUnion(m).cell.V1d),...
                   obj.meshUnion(m).cell.Npz, obj.meshUnion(m).cell.Nph, ...
                   obj.meshUnion(m).cell.V1d);
            end
        end
        
        function fphys = matLimitNew( obj, physClass, fphys )
            fphys{1}(:,:,physClass.varFieldIndex) = ...
                mxVertLimit3dNew(fphys{1}, physClass.varFieldIndex, obj.mesh.cell.Np, ...
                obj.mesh.K, obj.mesh.J, obj.mesh.cell.wq, obj.mesh.cell.Vq, physClass.fext3d{1},...
                physClass.gra, physClass.hcrit, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                int8(physClass.meshUnion.BottomBoundaryEdge.ftype), ...
                int8(physClass.meshUnion.SurfaceBoundaryEdge.ftype), obj.BoundaryEdge, ...
                obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, obj.Nv, obj.Nvc, obj.VToK,...
                physClass.meshUnion.mesh2d.cell.Fmask, physClass.meshUnion.mesh2d.cell.Np,...
                physClass.meshUnion.cell.Nz, physClass.meshUnion.mesh2d.cell.Nv);
        end
    end
    
end

