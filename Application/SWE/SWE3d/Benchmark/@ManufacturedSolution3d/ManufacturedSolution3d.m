classdef ManufacturedSolution3d < SWEBarotropic3d
    %MANUFACTUREDSOLUTION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 100;
        ChWidth = 100;
    end
    
    methods
        function obj = ManufacturedSolution3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3 10];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            %             obj.dt = 0.02;
            obj.Cf{1} = 0.005*ones(size(obj.mesh2d(1).x));
            
            obj.SurfBoundNewmannDate(:,:,1) = 1.5/1000 * ones(size(obj.SurfBoundNewmannDate(:,:,1)));%0.1
            %             Index =( all( obj.mesh2d.x - obj.ChLength/2  + 2*M > -1e-5 ));
            %             obj.WindTaux{1}(:,Index) = 0;
            %             Index =( all(obj.mesh2d.x + obj.ChLength/2 - 2*M < 1e-5 ));
            %             obj.WindTaux{1}(:,Index) = 0;
            %             obj.Cf{1} = 0.0025/1000;
        end
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            e = 0.001; w = 0.01; d = 0.1;
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                z = mesh3d.z; x = mesh3d.x; y = mesh3d.y;
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -(2-0.005*(mesh2d.x+mesh2d.y));
                
                fphys2d{m}(:,:,5) = ( mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                
                fphys2d{m}(:,:,6) = ( mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                %water depth
                fphys2d{m}(:,:,1) = e*(sin(w*(mesh2d.x+0)) + sin(w*(mesh2d.y+0))) - ...
                    fphys2d{m}(:, :, 4);
                
                fphys{m}(:,:,4) = mesh3d.Extend2dField(fphys2d{m}(:,:,1));
                %> Z in extended three dimensional fields
                fphys{m}(:,:,6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                %> '$\eta$' in extended three dimensional fields
                fphys{m}(:,:,7) = fphys{m}(:,:,4) + fphys{m}(:,:,6);
                Etax = mesh3d.rx .* (mesh3d.cell.Dr * fphys{m}(:,:,7)) + ...
                    mesh3d.sx .* (mesh3d.cell.Ds * fphys{m}(:,:,7));
                Etay = mesh3d.ry .* (mesh3d.cell.Dr * fphys{m}(:,:,7)) + ...
                    mesh3d.sy .* (mesh3d.cell.Ds * fphys{m}(:,:,7));
                
                fphys{m}(:,:,1) = -fphys{m}(:,:,4).*d.*z.*sin(w*(x+t));
                fphys{m}(:,:,2) = -fphys{m}(:,:,4).*d.*z.*sin(w*(y+t));
                
                fphys2d{m}(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys{m}(:, :, 1) );
                fphys2d{m}(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys{m}(:, :, 2) );
                
                fphys{m}(:,:,3) = -e*w.*cos(w.*(x+0)).*z - e*w.*cos(w.*(y+0)).*z + ...
                    fphys{m}(:,:,4).^2 .* (z.^2./2)*d*w.*cos(w.*(x+0)) - fphys{m}(:,:,4) .* d .*z .* sin(w.*(x+0)).*Etax + ...
                    fphys{m}(:,:,4).^2 .* (z.^2./2)*d*w.*cos(w.*(y+0)) - fphys{m}(:,:,4) .* d .*z .* sin(w.*(y+0)).*Etay;
                
            end
        end
        
        
        function [ option ] = setOption( obj, option )
            ftime = 50;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 100;
        end
        
        
        
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

% bctype = [ ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.ZeroGrad, ...
%     enumBoundaryCondition.ZeroGrad ];

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [ 0, obj.ChWidth ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end