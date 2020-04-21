classdef SWE3DVerticalExtendValue < SWE3DAbstractTest & ...
        SWEAbstract3d
    %SWE3DVERTICALEXTENDVALUE 此处显示有关此类的摘要
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
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    methods
        function obj = SWE3DVerticalExtendValue(N, Nz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz);
        end 
        
        function VerticalExtendValueTest(obj)
            Date = obj.mesh3d.Extend2dField( obj.mesh2d.x );
            obj.Assert(Date,obj.mesh3d.x);
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
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformTriMesh( N, ...
    [ 0, 10 ], [ 0, 10 ], 1, 4, bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, 4 );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

