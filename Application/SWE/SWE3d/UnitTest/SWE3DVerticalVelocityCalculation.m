classdef SWE3DVerticalVelocityCalculation < SWE3DAbstractTest & ...
        SWEBarotropic3d
    %SWE3DVERTICALINTEGRAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties

    end
    
    properties(Constant)
        hcrit = 0.001
    end
    
    methods
        function obj = SWE3DVerticalVelocityCalculation(N, Nz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz);
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
        end
        
        function VerticalIntegralTest(obj)
            
            Date = obj.mesh3d.VerticalColumnIntegralField( obj.fphys{1}(:, :, 1) );
            ExactDate =( 2 * 6 + 2 * obj.mesh2d.x )/6;
%             obj.Assert(ExactDate,Date);
            
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
            
            VerticalVelocity = obj.matEvaluateVerticalVelocity( obj.mesh3d, obj.fphys2d, obj.fphys );
%             obj.Assert(ExactVerticalVelocity,VerticalVelocity);
            
            VerticalVelocity = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, obj.fphys2d, obj.fphys );
            
            obj.Assert(ExactVerticalVelocity,VerticalVelocity);
            
            
 % For this case, vertical distribution of the first field changes abruptly  at position -2/6 and -4/6              
            index = find(all( abs(obj.mesh3d(1).z + 3/6)<= 1/6 + 1e-8 ));
            fphys{1}(:,index,1) = fphys{1}(:,index,1) + 5;
            fphys2d{1}(:,:,2) = obj.mesh3d.VerticalColumnIntegralField( obj.fphys{1}(:, :, 1) );
            fphys2d{1}(:,:,3) = obj.mesh3d.VerticalColumnIntegralField( obj.fphys{1}(:, :, 2) );
            VerticalVelocity = obj.matEvaluateVerticalVelocity( obj.mesh3d, obj.fphys2d, obj.fphys );
            obj.Assert(ExactVerticalVelocity,VerticalVelocity);            
        end
        
    end
    methods(Access = protected)
        
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d{1} = zeros(size(obj.mesh2d.x,1),size(obj.mesh2d.x,2),obj.Nfield2d);
            fphys{1}   = zeros(size(obj.mesh3d.x,1),size(obj.mesh3d.x,2),obj.Nfield);
            fphys{1}(:,:,1) = 2 * ones(size(obj.mesh3d.x));
 % For this case, vertical distribution of the first field is continuously distributed           
            index = (obj.mesh3d(1).z==-0.5);
            fphys{1}(index) = fphys{1}(index) + 2 * obj.mesh3d.x(index);
            
            fphys2d{1}(:,:,1) = 1 * ones(size(obj.mesh2d.x));
            fphys2d{1}(:, :, 4) = -1 * fphys2d{1}(:,:,1);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 21600;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;                  
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantEddyViscosityValue') = 0.01;
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
    [ 0, 10 ], [ 0, 10 ], 2, 2, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, 6 );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, 6 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, 6 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

