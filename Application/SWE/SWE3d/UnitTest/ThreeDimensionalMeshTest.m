classdef ThreeDimensionalMeshTest < SWE3DAbstractTest & ...
        SWEAbstract3d
    %THREEDIMENSIONALMESHTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Nfield2d = 1
        Nvar2d = 1
        varFieldIndex2d = 1
        outputFieldOrder2d = []
        outputFieldOrder =  1
        Nfield = 1
        Nvar = 1
        varFieldIndex = 1
        N
        Nz
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    methods
        function obj = ThreeDimensionalMeshTest(N, Nz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz);
            obj.N = N;
            obj.Nz = Nz;
        end
        
        function QuadMeshMatrixTest(obj)
            % for this test,the mesh is carefully designed
            % used to verify the equality between mass matrix of standard cell and that of edge
            
            obj.Assert(obj.mesh2d.cell.M, obj.mesh3d.InnerEdge.M);
            obj.Assert(obj.mesh2d.cell.M, obj.mesh3d.BoundaryEdge.M);
            obj.Assert(obj.mesh2d.cell.M, obj.mesh3d.BottomBoundaryEdge.M);
            obj.Assert(obj.mesh2d.cell.M, obj.mesh3d.SurfaceBoundaryEdge.M);
            % used to verify the equality between Jacobian of two-dimensional mesh and that of edge             
            obj.Assert(obj.mesh2d.J(:,1), obj.mesh3d.InnerEdge.Js(:,1));
            obj.Assert(obj.mesh2d.J(:,1), obj.mesh3d.BoundaryEdge.Js(:,1));
            obj.Assert(obj.mesh2d.J(:,1), obj.mesh3d.BottomBoundaryEdge.Js(:,1));
            obj.Assert(obj.mesh2d.J(:,1), obj.mesh3d.SurfaceBoundaryEdge.Js(:,1));
        end
        
        function TriMeshMatrixTest(obj)
            bctype = [ ...
                enumBoundaryCondition.ZeroGrad, ...
                enumBoundaryCondition.ZeroGrad, ...
                enumBoundaryCondition.ZeroGrad, ...
                enumBoundaryCondition.ZeroGrad ];
            
            mesh2d = makeUniformTriMesh( obj.N, ...
                [ -1, 0 ], [ -1, 0 ], 5, 5, bctype);
            
            cell = StdPrismTri( obj.N, obj.Nz );
            zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
            mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, 5 );
            mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
            mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
            mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
            mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
            mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
            % used to verify the equality between mass matrix of standard cell and that of edge             
            obj.Assert(obj.mesh2d.cell.M, mesh3d.InnerEdge.M);
            obj.Assert(obj.mesh2d.cell.M, mesh3d.BoundaryEdge.M);       
        end
    end
    
    methods( Hidden )
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            %an empty function
        end
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz )

bctype = [ ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -1, 0 ], [ -1, 0 ], 5, 5, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, 5 );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

