classdef TideInSemiClosedSlopeChannel3d< SWEBarotropic3d
    %TIDEINSEMICLOSEDCHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties( Constant )
        
        SMSFile = [pwd,'\Application\SWE\SWE3d\Benchmark\@TideInSemiClosedSlopeChannel3d\Data\fort.14'];
        
    end
    
    properties(Constant)
        amplitude = 0.5
        H0 = 6
        hcrit = 0.01;
        T = 12*3600
        %         omega = 2*pi/12/3600
        ChLength = 3400
        ChWidth = 200
        startTime = 0
        finalTime = 5*3600
    end
    
    properties
        Length
        k
    end
    
    methods
        function obj = TideInSemiClosedSlopeChannel3d( N, Nz, M, Mz )
            % setup mesh domain
%             [ mesh2d, mesh3d ] = obj.makeChannelMesh( N, Nz, Mz );
            [ mesh2d, mesh3d ] = obj.makeChannelMesh( N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 4];
            obj.outputFieldOrder3d = [ 1 2 3 10];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( mesh2d, mesh3d );
            %> time interval
            obj.Cf{1} = 0;
        end
        
    end
    
    methods(Access = protected)
        
        function matUpdateExternalField( obj, time, ~, ~ )
            
            obj.fext2d{1}( :, :, 3 ) = obj.amplitude * cos(2*pi/obj.T*time + pi/2 )+obj.H0;
            obj.fext3d{1}( :, :, 3 ) = obj.amplitude * cos(2*pi/obj.T*time + pi/2 )+obj.H0;
            
        end
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
%             fphys2d = cell( obj.Nmesh, 1 );
%             fphys = cell( obj.Nmesh, 1 );
%             fphys2d{1} = zeros( obj.mesh2d(1).cell.Np, obj.mesh2d(1).K, obj.Nfield2d );
%             fphys2d = obj.matInterpolateTopography(fphys2d);
%             fphys2d{1}(:,:,1) = -1*fphys2d{1}(:,:,4);
%             fphys{1} = zeros( obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nfield );

            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            fphys2d{1} = zeros( obj.mesh2d(1).cell.Np, obj.mesh2d(1).K, obj.Nfield2d );
            fphys2d{1}(:,:,4) = -1 + (-6 + 1)/obj.ChLength * obj.mesh2d(1).x;
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
            outputIntervalNum = 1500;
            option('startTime') = obj.startTime;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = (ftime - obj.startTime)/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 230;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('EosType') = enumEOSType.Linear;
            option('ConstantVerticalEddyViscosityValue') = 0.1;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 70;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.None;
            %             option('f0 for beta coriolis solver') = 0;
            option('beta for beta coriolis solver') = 0.0;
            option('PhysicalSurfaceRoughnessLength') = 0.0001;
            option('PhysicalBottomRoughnessLength') = 0.0001;
        end
        
        function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )
            
            bctype = [ ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.ClampedDepth ];
            
            mesh2d = makeUniformQuadMesh(N, [0, obj.ChLength],...
                [-obj.ChWidth/2, obj.ChWidth/2], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);
            
            cell = StdPrismQuad( N, Nz );
            zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
            mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
            mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
            mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
            mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
            
        end
        
        %         function  [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, Mz)
        %             [ mesh2d ] = makeSMSFileUMeshUnion2d( N, obj.SMSFile );
        %             cell = StdPrismTri( N, Nz );
        %             zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
        %             mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
        %             mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
        %             mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
        %             mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
        %             mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
        %             mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
        %         end
        
    end
    
end


