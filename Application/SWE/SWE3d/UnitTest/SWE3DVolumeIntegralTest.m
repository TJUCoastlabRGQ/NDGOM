classdef SWE3DVolumeIntegralTest < SWE3DAbstractTest & ...
        SWEBarotropic3d
    %SWE3DFLUXCONSISTENCYTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        hcrit = 1
    end    
    
    methods
        function obj = SWE3DVolumeIntegralTest(N, Nz, Mz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz, Mz);
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            obj.SurfBoundNewmannDate(:,:,1) = zeros(size(obj.mesh2d(1).x));
        end
        
        function VolumeIntegralTest(obj)
            % Results indicate these two are not equivalent
            obj.fphys{1}(:,:,1) = 0.5 * rand(size(ones(size(obj.mesh3d.x))));
            obj.fphys{1}(:,:,2) = 0.5 * rand(size(ones(size(obj.mesh3d.x))));
            
            obj.fphys2d{1}(:,:,2) = obj.mesh3d.VerticalColumnIntegralField( obj.fphys{1}(:, :, 1) );
            obj.fphys2d{1}(:,:,3) = obj.mesh3d.VerticalColumnIntegralField( obj.fphys{1}(:, :, 2) );
            
            Derive3d = obj.mesh3d.VerticalColumnIntegralField( obj.mesh3d.rx .* (obj.mesh3d.cell.Dr * (obj.fphys{1}(:,:,1))) + obj.mesh3d.sx .* (obj.mesh3d.cell.Ds * (obj.fphys{1}(:,:,1))) + ...
               obj.mesh3d.ry .* (obj.mesh3d.cell.Dr * (obj.fphys{1}(:,:,2))) + obj.mesh3d.sy .* (obj.mesh3d.cell.Ds * (obj.fphys{1}(:,:,2))) );
           
           Derive2d = obj.mesh2d.rx .* (obj.mesh2d.cell.Dr * obj.fphys2d{1}(:,:,2)) + obj.mesh2d.sx .* (obj.mesh2d.cell.Ds * obj.fphys2d{1}(:,:,2)) + ...
               obj.mesh2d.ry .* (obj.mesh2d.cell.Dr * obj.fphys2d{1}(:,:,3)) + obj.mesh2d.sy .* (obj.mesh2d.cell.Ds * obj.fphys2d{1}(:,:,3));
            
            obj.Assert(Derive3d, Derive2d);
            
        end        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
            end            
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

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, 10 ], [ 0, 10 ], 5, 5, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end