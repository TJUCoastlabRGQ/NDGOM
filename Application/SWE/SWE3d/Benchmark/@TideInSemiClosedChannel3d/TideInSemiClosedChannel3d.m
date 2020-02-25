classdef TideInSemiClosedChannel3d < SWEBarotropic3d
    %TIDEINSEMICLOSEDCHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        amplitude = 0.25
        H0 = 30
        hcrit = 0.01;
        T = 12*3600
%         omega = 2*pi/12/3600
        ChLength = 24000
        ChWidth = 2000
        
    end
    
    properties
        Length
        k
        dt
    end
    
    methods
        function obj = TideInSemiClosedChannel3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3];
            obj.outputFieldOrder = [ 1 2 3 4];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            obj.dt = 0.02;
            obj.Cf{1} = 0;
            obj.WindTaux{1} = zeros(size(obj.mesh2d(1).x));
            obj.WindTauy{1} = zeros(size(obj.mesh2d(1).y));            
%             obj.Cf{1} = 0.0025/1000;
        end
        
    end
    
    methods(Access = protected)

        function WaveCharacterEstimate(obj)
            f = @(L) L - obj.gra*(obj.T)^2/(2*pi)*tanh(2*pi/L*obj.Depth);
            obj.Length = fzero(f,[0.01 20000000]);
        end
        
        function matUpdateExternalField( obj, time, ~ )
            
            obj.fext{1}( :, :, 1 ) = obj.amplitude * sin(2*pi/obj.T*time )+obj.Depth;
            
        end
        
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;                
                %water depth
                fphys2d{m}(:,:,1) = obj.H0;
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
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('EddyViscosityType') = enumEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantEddyViscosityValue') = 0;
        end
        
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedDepth ];

mesh2d = makeUniformTriMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

