classdef NdgGaussQuadWeakFormSolver3d < NdgGaussQuadWeakFormSolver
    %NDGGAUSSQUADWEAKFORMSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %> SBLIFT, values of the basis functions at each edge quadrature points of the surface boundary edges
        %> \f$$ [\mathcal{LIFT}_t]_{i,ej} = \varphi_i(\mathbf{r}_{ej}), \f$$
        %> where \f$ \mathbf{r}_{ej} \f$ is the location of quadrature points on edges.
        %> Likewise, BBLIFT stands for that at each edge quadrature points of the bottom boundary edges,
        %> and BELIFT stands for that at each edge quadrature points of the bottom edges,
        SBLIFT
        BOTLIFT
        BBLIFT
        %> outward normal vector at the local surface boundary edge
        SBnx, SBny, SBnz
        %> outward normal vector at the bottom edge
        BOTnx, BOTny, BOTnz
        %> outward normal vector at the bottom boundary edge
        BBnx, BBny, BBnz
        %> the element in the vector equals to
        %> \f$ w_{ei} \cdot Js_{ei} \f$
        %> where \f$ w_{ei} \f$ and \f$ Js_{ei} \f$ are the surface integral weights and Jacobian
        %> at each quadrature points on local surface boundary edges.
        SBwJs
        %> Likewise, BBwJs are the surface integral weights and Jacobian at each quadrature points on bottom boundary edges,
        %> and BEwJs are the surface integral weights and Jacobian at each quadrature points on bottom edges,
        BOTwJs
        BBwJs
        %> project nodal values to face quadrature points
        %         Vfq
        %> project face values to face quadrature points at the surface boundary edges
        SBFVfq
        %> project face values to face quadrature points at the bottom edges
        BOTFVfq
        %> project face values to face quadrature points at the bottom boundary edges
        BBFVfq
    end
    
    methods
        function obj = NdgGaussQuadWeakFormSolver3d( phys, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver(phys, meshUnion );
            for m = 1:phys.Nmesh
                mesh = meshUnion( m );
                %                 stdcell = mesh.cell;
                [ obj.SBFVfq{m}, obj.BOTFVfq{m}, obj.BBFVfq{m}, ...
                    obj.SBwJs{m}, obj.BOTwJs{m}, obj.BBwJs{m}] = obj.assemble3dFacialVandMatrixFaceQuadrature( mesh );
                [ obj.SBnx{m}, obj.SBny{m}, obj.SBnz{m},...
                    obj.BOTnx{m}, obj.BOTny{m}, obj.BOTnz{m},...
                    obj.BBnx{m}, obj.BBny{m}, obj.BBnz{m} ] = obj.assemble3dNormalVector( mesh, obj.SBFVfq{m}, obj.BOTFVfq{m}, obj.BBFVfq{m} );
                [ obj.SBLIFT{m}, obj.BOTLIFT{m}, obj.BBLIFT{m} ] = obj.assemble3dBoundaryLiftMatrix( mesh );
            end
        end
    end
    
    methods( Static )
        function [ SBFVfq, BOTFVfq, BBFVfq, SBwJs, BOTwJs, BBwJs] = assemble3dFacialVandMatrixFaceQuadrature( mesh )
            fcell = getStdCell( mesh.cell.N, mesh.SurfaceBoundaryEdge.type );
            SBFVfq = fcell.nodal_func(fcell.rq, fcell.sq, fcell.tq);
            BOTFVfq = SBFVfq; BBFVfq = SBFVfq;
            SBwJs = bsxfun(@times, fcell.wq, SBFVfq * mesh.SurfaceBoundaryEdge.Js);
            BOTwJs = SBwJs; BBwJs = SBwJs;
        end
        
        function  [ SBnx, SBny, SBnz, BOTnx, BOTny, BOTnz, BBnx, BBny, BBnz ] = assemble3dNormalVector( mesh, SBFVfq, BOTFVfq, BBFVfq )
            SBnx = SBFVfq * mesh.SurfaceBoundaryEdge.nx; BOTnx = BOTFVfq * mesh.BottomEdge.nx; BBnx = BBFVfq * mesh.BottomBoundaryEdge.nx;
            SBny = SBFVfq * mesh.SurfaceBoundaryEdge.ny; BOTny = BOTFVfq * mesh.BottomEdge.ny; BBny = BBFVfq * mesh.BottomBoundaryEdge.ny;
            SBnz = SBFVfq * mesh.SurfaceBoundaryEdge.nz; BOTnz = BOTFVfq * mesh.BottomEdge.nz; BBnz = BBFVfq * mesh.BottomBoundaryEdge.nz;
        end
        
        function [ SBLIFT, BOTLIFT, BBLIFT ] = assemble3dBoundaryLiftMatrix( mesh )
            fcell = getStdCell( mesh.cell.N, mesh.SurfaceBoundaryEdge.type );
            SBLIFT = ( fcell.nodal_func(fcell.rq, fcell.sq, fcell.tq) )';
            BOTLIFT = SBLIFT; BBLIFT = SBLIFT;
        end
        
    end
    
end

