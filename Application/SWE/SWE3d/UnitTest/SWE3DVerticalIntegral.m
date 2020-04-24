classdef SWE3DVerticalIntegral < SWE3DAbstractTest & ...
        SWEAbstract3d
    %SWE3DVERTICALINTEGRAL 此处显示有关此类的摘要
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
        function obj = SWE3DVerticalIntegral(N, Nz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz);
        end
        
        function VerticalIntegralTest(obj)
            fphys = cell(1);
            fphys{1} = zeros(size(obj.mesh3d.x));
            index = (obj.mesh3d(1).z>=-0.5);
            fphys{1}(index) = -12/5*obj.mesh3d.z(index);
            fphys{1}(~index) = 12/5*obj.mesh3d.z(~index)+2.4;
            Date = obj.mesh3d.VerticalColumnIntegralField( fphys{1}(:, :, 1) );
            ExactDate = 0.6 * ones(size(obj.mesh2d.x));
            obj.Assert(ExactDate,Date);
            
            index = (obj.mesh3d.x ==0 & obj.mesh3d.y == 0);
            fphys{1}(index) = fphys{1}(index) + 2;
            Date = obj.mesh3d.VerticalColumnIntegralField( fphys{1}(:, :, 1) );
            ExactDate(1) = ExactDate(1) + 2;
            obj.Assert(ExactDate,Date);
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

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, 10 ], [ 0, 10 ], 1, 4, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, 4 );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

