classdef SWE3DVerticalVelocityCalculation < SWE3DAbstractTest & ...
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
        function obj = SWE3DVerticalVelocityCalculation(N, Nz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz);
        end
        
        function VerticalIntegralTest(obj)
            fphys = cell(1);
            fphys{1}(:,:,1) = 2 * ones(size(obj.mesh3d.x));
            
            index = (obj.mesh3d(1).z==-0.5);
            fphys{1}(index) = fphys{1}(index) + 2 * obj.mesh3d.x(index);
            
            Date = obj.mesh3d.VerticalColumnIntegralField( fphys{1}(:, :, 1) );
            ExactDate =( 2 * 6 + 2 * obj.mesh2d.x )/6;
            obj.Assert(ExactDate,Date);
            
            ExactVerticalVelocity = zeros(size(obj.mesh3d.x));
            index = (abs(obj.mesh3d(1).z + 5/6) < 1e-8 );
            ExactVerticalVelocity(index) = 1/18;
            index = (abs(obj.mesh3d(1).z + 4/6) < 1e-8 );
            ExactVerticalVelocity(index) = 1/9;    
            index = (abs(obj.mesh3d(1).z + 3/6) < 1e-8 );
            ExactVerticalVelocity(index) = 0;  
            index = (abs(obj.mesh3d(1).z + 2/6) < 1e-8 );
            ExactVerticalVelocity(index) = -1/9;  
            index = (abs(obj.mesh3d(1).z + 1/6) < 1e-8 );
            ExactVerticalVelocity(index) = -1/18;
            index = (abs(obj.mesh3d(1).z)  < 1e-8 );
            ExactVerticalVelocity(index) = 0; 
            
            fphys2d = cell(1);
            fphys{1}(:,:,2) = zeros(size(fphys{1}(:,:,1)));
            fphys2d{1}(:,:,1) = zeros(size(obj.mesh2d.x));
            fphys2d{1}(:,:,2) = obj.mesh3d.VerticalColumnIntegralField( fphys{1}(:, :, 1) );
            fphys2d{1}(:,:,3) = obj.mesh3d.VerticalColumnIntegralField( fphys{1}(:, :, 2) );
            VerticalVelocity = obj.matEvaluateVerticalVelocity( obj.mesh3d, fphys2d{1}, fphys{1} );
            obj.Assert(ExactVerticalVelocity,VerticalVelocity);
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

mesh2d = makeUniformTriMesh( N, ...
    [ 0, 10 ], [ 0, 10 ], 1, 6, bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, 6 );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

