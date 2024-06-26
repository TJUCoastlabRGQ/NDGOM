classdef SWE3DFluxConsistencyTest < SWE3DAbstractTest & ...
        SWEBarotropic3d
    %SWE3DFLUXCONSISTENCYTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
%         Nfield2d = 4
%         Nvar2d = 1
%         varFieldIndex2d = 1
%         outputFieldOrder2d = []
%         outputFieldOrder =  1
%         Nfield = 9
%         Nvar = 1
%         varFieldIndex = 1
    end
    
    properties(Constant)
        hcrit = 1
    end    
    
    methods
        function obj = SWE3DFluxConsistencyTest(N, Nz, Mz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz, Mz);
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            obj.SurfBoundNewmannDate(:,:,1) = zeros(size(obj.mesh2d(1).x));
        end
        
        function FluxConsistencyTest(obj)
            % Results indicate these two are not equivalent
            obj.fphys{1}(:,:,1) = 0.5 * rand(size(ones(size(obj.mesh3d.x))));
%             obj.fphys{1}(:,:,1) = 0.5 * ones(size(obj.mesh3d.x));
            index = (obj.mesh3d(1).z==-0.5);
            obj.fphys{1}(index) = 0.75;
            
            obj.fphys{1}(:,:,2) = 0.5 * rand(size(ones(size(obj.mesh3d.x))));
%             obj.fphys{1}(:,:,2) = 0.5 * ones(size(obj.mesh3d.x));
%             obj.fphys{1}(5 + numel(obj.mesh3d.x)) = 0.9;
            
            obj.fphys2d{1}(:,:,1) = ( 10 + rand(1) ) * ones(size(obj.mesh2d.x)); 
%             obj.fphys2d{1}(:,:,1) = 10 * ones(size(obj.mesh2d.x)); 
            obj.fphys{1}(: , :, 4) = obj.mesh3d.Extend2dField( obj.fphys2d{1}(:, :, 1) );
            
            obj.fphys2d{1}(:,:,2) = obj.mesh3d.VerticalColumnIntegralField( obj.fphys{1}(:, :, 1) );
            obj.fphys2d{1}(:,:,3) = obj.mesh3d.VerticalColumnIntegralField( obj.fphys{1}(:, :, 2) );
            
            InnerEdge2d = obj.mesh2d.InnerEdge;
            [ fm2d, fp2d ] = InnerEdge2d.matEvaluateSurfValue( obj.fphys2d );
            FluxS2d = obj.matEvaluateSurfNumFlux(obj.mesh2d, InnerEdge2d.nx, InnerEdge2d.ny, fm2d, fp2d, InnerEdge2d);
            
            InnerEdge3d = obj.mesh3d.InnerEdge;
            [ fm3d, fp3d ] = InnerEdge3d.matEvaluateSurfValue( obj.fphys );
            FluxS3d = obj.matEvaluateSurfNumFlux(obj.mesh3d, InnerEdge3d.nx, InnerEdge3d.ny, fm3d(:,:,[4,1,2]), fp3d(:,:,[4,1,2]), InnerEdge3d);
            
            obj.Assert(InnerEdge3d.VerticalColumnIntegralField(FluxS3d(:,:,1)),FluxS2d(:,:,1));
            
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
    [ 0, 10 ], [ 0, 10 ], 3, 3, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end