classdef WellBalanceTest < SWEBarotropic3d
    %WINDDRIVENFLOW 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )

        startTime = 0;
        finalTime = 60;
        hcrit = 0.001;
        
    end
    
    methods
        function obj = WellBalanceTest( N, Nz, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeMesh( obj, N, Nz, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3 10];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            obj.Cf{1} = 0.005*ones(size(obj.mesh2d(1).x));
            
            obj.SurfBoundNewmannDate(:,:,1) = 0.0/1000 * ones(size(obj.SurfBoundNewmannDate(:,:,1)));%0.1
        end
        
        matExplicitRK222(obj);
        
        
        EntropyAndEnergyCalculation(obj);
        
        AnalysisResult2d( obj );
        AnalysisResult3d( obj );
        plot_uv2d( obj );
        plot_h2d( obj );
        PostProcess_level( obj );
        PostProcess_speed( obj );
        PostProcess_direction( obj );
    end
    
    methods ( Access = protected )
        
        %> set initial function
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            fphys2d{1} = zeros( obj.mesh2d(1).cell.Np, obj.mesh2d(1).K, obj.Nfield2d );
            fphys2d = obj.matInterpolateTopography(fphys2d);
            fphys2d{1}(:,:,1) = 0.3 -1*fphys2d{1}(:,:,4);
            
%             Index = find(obj.meshUnion.mesh2d.x >500005 & obj.meshUnion.mesh2d.y < 4306120);
%             fphys2d{1}(3*obj.meshUnion.mesh2d.K*obj.meshUnion.mesh2d.cell.Np + Index) = -28.0;

            
%             Index = find(abs(obj.meshUnion.mesh2d.x - 498357)<1 & abs(obj.meshUnion.mesh2d.y - 4311231)<1);
%             fphys2d{1}(3*obj.meshUnion.mesh2d.K*obj.meshUnion.mesh2d.cell.Np + Index) = -18.0;
            
            fphys{1} = zeros( obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nfield );
        end
        
         function fphys2d = matInterpolateTopography( obj, fphys2d )
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            for m = 1:obj.Nmesh
                mesh = obj.mesh2d(m);
                tpfile = [path, '/private/Bathmetry.txt'];
                fid1 = fopen(tpfile, 'r');
                % READ FIRST LINE
                fgetl(fid1);
                % READ POINT
                temp = fscanf(fid1,'%d',2); Nv = temp(2);
                data = fscanf(fid1,'%d %f %f %f\n',[4,Nv]);
                fclose(fid1);
                interp = scatteredInterpolant( ...
                    data(2,:)', data(3,:)', data(4,:)', 'linear');
                
                fphys2d{m}(:,:,4) = interp(mesh.x, mesh.y);
            end
        end% func       
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
% %             obj.BotBoundNewmannDate(:,:,1)  = obj. 
%            VCV = obj.meshUnion(1).cell.VCV;
%            Nz = obj.meshUnion(1).Nz;
%            Hu = VCV * fphys{1}(:,Nz:Nz:end,1);
%            Hv = VCV * fphys{1}(:,Nz:Nz:end,2);
%            H  = VCV * fphys{1}(:,Nz:Nz:end,4);
%            obj.BotBoundNewmannDate(:,:,1) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
%                (Hv./H).^2 ) .* ( Hu./H ) * (-1);
%            obj.BotBoundNewmannDate(:,:,2) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
%                (Hv./H).^2 ) .* ( Hv./H ) * (-1);           
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 10800;
            outputIntervalNum = 3600;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.01;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 0.1;
        end
        
    end
    
end

function [mesh2d, mesh3d] = makeMesh( obj, N, Nz, Mz )

% bctype = [ ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.ZeroGrad, ...
%     enumBoundaryCondition.ZeroGrad ];


mesh2d = makeSMSFileUMeshUnion2d( N, [pwd,'\Application\SWE\SWE3d\Benchmark\@WellBalanceTest\private\fort_z0.14'] );

cell = StdPrismTri( N, Nz );
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

