classdef BohaiEstuary < SWEBarotropic3d
    %BOHAIESTUARY 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties( Constant )
        
        hcrit = 0.1;
        
        SMSFile = [pwd,'\Application\SWE\SWE3d\Benchmark\@BohaiEstuary\Data\fort.14'];
        
        tidalFile = [pwd,'\Application\SWE\SWE3d\Benchmark\@BohaiEstuary\Data\OutputBoundary.xlsx'];
        
        GotmFile = [pwd,'\Application\SWE\SWE3d\Benchmark\@BohaiEstuary\Data\gotmturb.nml'];
    end
    
    properties
        finalTime = 1435*600 + 38*3600;
        
        tideinterval = 600;
        
        OBEdgeIndex
        
        Tide
    end
    
    properties
        T0
        
        S0
        
        alphaT
        
        betaS
    end    
    
    methods
        function obj = BohaiEstuary(N, Nz, Mz)
            %BOHAIESTUARY 构造此类的实例
            %   此处显示详细说明
            [ mesh2d, mesh3d ] = obj.makeChannelMesh( N, Nz, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 4 6 ];
            obj.T0 = 0;
            obj.S0 = 0;
            obj.alphaT = 0;
            obj.betaS = 0;
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.ReadTideElevation;
        end
    end
    
    methods( Access = protected)
        
        ReadTideElevation( obj );
        
        function boundaryInfo = initBoundaryInformation( obj )
            boundaryInfo = struct('time', {}, 'extElevation', {});
            
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            obcfile = [ path, '/mesh/Benchmark_2_input.txt'];
            fp = fopen(obcfile);
            fgetl(fp); % pass the rest of first line
            data = fscanf(fp, '%f %f', [2, inf]);
            fclose(fp);
            
            boundaryInfo(1).time = data(1, :);
            boundaryInfo(1).extElevation = data(2, :);
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                edge = obj.meshUnion(m).BoundaryEdge;
                nodeid = bsxfun( @plus, edge.FToN1, (edge.FToE(1, :) - 1) .* mesh.cell.Np);
                obj.fext{m}( :, :, 4 ) = obj.fphys{m}( nodeid + mesh.K * mesh.cell.Np * 3 );
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
            delta = obj.tideinterval;
            for m1 = 1:obj.Nmesh
                %     mesh = obj.meshUnion( m1 );
                
                % time step
                s1 = floor( time/delta ) + 1;
                s2 = s1 + 1;
                alpha1 = ( delta * (s2 - 1) - time ) / delta;
                alpha2 = ( time - delta * (s1 - 1) ) / delta;
                
                ind = obj.OBEdgeIndex{m1};
                fnT = obj.Tide{m1}(:, :, s1) .* alpha1 + obj.Tide{m1}(:, :, s2) .* alpha2;
                obj.fext2d{m1}(:,ind,1) = max( fnT - obj.fext2d{m1}(:,ind,4), 0 );
                obj.fext3d{m1}(:,:,3) = obj.meshUnion.BoundaryEdge.matVerticalRepmatFacialValue(obj.fext2d{m1}(:,:,1));
            end
        end
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            fphys2d{1} = zeros( obj.mesh2d(1).cell.Np, obj.mesh2d(1).K, obj.Nfield2d );
            fphys2d = obj.matInterpolateTopography(fphys2d);
            fphys2d{1}(:,:,1) = -1*fphys2d{1}(:,:,4);
            fphys{1} = zeros( obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nfield );
        end
        
        function fphys2d = matInterpolateTopography( obj, fphys2d )
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            for m = 1:obj.Nmesh
                mesh = obj.mesh2d(m);
                tpfile = [path, '/Data/Bathmetry.txt'];
                fid1 = fopen(tpfile, 'r');
                % READ FIRST LINE
                fgetl(fid1);
                % READ POINT
                temp = fscanf(fid1,'%d',2); Nv = temp(2);
                data = fscanf(fid1,'%d %f %f %f\n',[4,Nv]);
                fclose(fid1);
                interp = scatteredInterpolant( ...
                    data(2,:)', data(3,:)', -1*data(4,:)', 'linear');
                
                fphys2d{m}(:,:,4) = interp(mesh.x, mesh.y);
            end
        end% func
        
        function [ option ] = setOption( obj, option )
            ftime = obj.finalTime;
            outputIntervalNum = 440;
            option('startTime') = 1435*600;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 400;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('limiterType') = enumLimiter.Vert;
            option('EosType') = enumEOSType.Linear;
            option('ConstantVerticalEddyViscosityValue') = 0.1;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.Constant;
            option('ConstantHorizontalEddyViscosityValue') = 70;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.Beta;
            option('f0 for beta coriolis solver') = 2*sin((38 + 54/60)/180*pi)*7.29*10^(-5);
            %             option('f0 for beta coriolis solver') = 0;
            option('beta for beta coriolis solver') = 0.0;
            option('PhysicalSurfaceRoughnessLength') = 0.003;
            option('PhysicalBottomRoughnessLength') = 0.003;
        end
        
        
        function  [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, Mz)
            [ mesh2d ] = makeSMSFileUMeshUnion2d( N, obj.SMSFile );
            cell = StdPrismTri( N, Nz );
            zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
            mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
            mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
            mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
            mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
        end
    end
end

