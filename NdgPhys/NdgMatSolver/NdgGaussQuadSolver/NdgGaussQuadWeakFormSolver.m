classdef NdgGaussQuadWeakFormSolver < NdgGaussQuadStrongFormSolver
    
    methods
        function obj = NdgGaussQuadWeakFormSolver( phys )
            obj = obj@NdgGaussQuadStrongFormSolver( phys );

            for m = 1:phys.Nmesh
                mesh = phys.meshUnion( m );
                [ obj.Vq{m} ] = mesh.cell.Vq;
                obj.Dr{m} = obj.Dr{m}';
                obj.Ds{m} = obj.Ds{m}';
                obj.Dt{m} = obj.Dt{m}';
            end
        end% func
    end
    
    methods ( Access = protected )
        %> Function used to interpolate the value from the interpolation point to the Volume Gauss quadrature point
        function InterVolumefphys = matInterpolateToVolumeGaussQuadraturePoint(obj, Vq, fphys)
            InterVolumefphys = bsxfun(@times, Vq, fphys);
        end
        
        function [ InterFaceLocalfphys, InterFaceAdjacentfphys ] = matInterpolateToFaceGaussQuadraturePoint(obj, edge, Vfq, fphys)
            for i = 1:size(fphys, 3)
                InterFaceLocalfphys = bsxfun(@times, Vfq, fphys( edge.GFToN1 + (i-1) * numel(fphys(:,:,i))) );
                InterFaceAdjacentfphys = bsxfun(@times, Vfq, fphys( edge.GFToN2 + (i-1) * numel(fphys(:,:,i))) );
            end
        end
    end
    
end