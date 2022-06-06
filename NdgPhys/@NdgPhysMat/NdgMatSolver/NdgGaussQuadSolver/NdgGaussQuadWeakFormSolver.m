classdef NdgGaussQuadWeakFormSolver < NdgGaussQuadStrongFormSolver
    
    methods
        function obj = NdgGaussQuadWeakFormSolver( phys, meshUnion )
            obj = obj@NdgGaussQuadStrongFormSolver( phys, meshUnion );
            
            for m = 1:phys.Nmesh
                mesh = meshUnion( m );
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
%                 InterVolumefphys= bsxfun(@times, Vq, fphys);
               InterVolumefphys = zeros(size(Vq,1), size(fphys,2), size(fphys, 3));
               for i = 1:size( fphys, 3 )
%                    tempfphys = fphys(:,:,i);
%                    InterVolumefphys(:,:,i)= bsxfun(@times, Vq, fphys(:,:,i));
                   InterVolumefphys(:,:,i)=  Vq* fphys(:,:,i);
               end
                
%               InterVolumefphys =  bsxfun(@times, Vq, permute( permute( fphys, [1,3,2] ), [2,1,3] ) );
        end
        
        function [ GaussInterFaceLocalfphys, GaussInterFaceAdjacentfphys ] = matInterpolateToFaceGaussQuadraturePoint(obj, edge, Vfq, InterFaceLocalfphys, InterFaceAdjacentfphys )
            GaussInterFaceLocalfphys = zeros(size(Vfq,1), size(InterFaceLocalfphys,2), size(InterFaceLocalfphys, 3));
            GaussInterFaceAdjacentfphys = zeros(size(Vfq,1), size(InterFaceLocalfphys,2), size(InterFaceLocalfphys, 3));
            for i = 1:size(InterFaceLocalfphys, 3)
%                 InterFaceLocalfphys = bsxfun(@times, Vfq, InterFaceLocalfphys(:,:,i) );
%                 InterFaceAdjacentfphys = bsxfun(@times, Vfq, InterFaceAdjacentfphys(:,:,i) );
                GaussInterFaceLocalfphys(:,:,i) =  Vfq * InterFaceLocalfphys(:,:,i);
                GaussInterFaceAdjacentfphys(:,:,i) =  Vfq * InterFaceAdjacentfphys(:,:,i);
            end
        end
        
        function RHS = matAssembleIntoRHS( obj, edge, FRHS, RHS)
            for i = 1:numel(FRHS)
                RHS(edge.GFToN1(i)) = RHS(edge.GFToN1(i)) + FRHS(i);
                RHS(edge.GFToN2(i)) = RHS(edge.GFToN2(i)) - FRHS(i);
            end
        end
 
        
        function RHS = matAssembleBoundaryAndSourceTermIntoRHS( obj, edge, FRHS, RHS)
            for i = 1:numel(FRHS)
                RHS(edge.GFToN1(i)) = RHS(edge.GFToN1(i)) + FRHS(i);
            end
        end
    end
    
end
