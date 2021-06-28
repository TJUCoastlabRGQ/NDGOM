classdef StandingWaveInAClosedChannel < SWEBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        hcrit = 0.01;
        %         ChLength = 100;
        ChLength = 20;
        %> channel width
        ChWidth = 2;
        %> channel depth
        H0 = 10;
        %> x range
        %> start time
        startTime = 0;
        %         %> final time
        finalTime = 10;
        % to be corrected
        GotmFile = fullfile('D:\PhdResearch\Application\SWE\SWE3d\Benchmark\@StandingWaveInAClosedChannel','\gotmturb.nml');
    end
    
    properties
        SurfaceBoundaryEdgeType = 'Dirichlet'
    end
    
    properties
        dt
        miu0
        Lambda = 20;
        A = 0.1;
    end
    
    methods
        function obj = StandingWaveInAClosedChannel( N, Nz, M, Mz )
            % setup mesh domain
            [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [1 2 3 11];
            obj.Nfield = 11;
            obj.Nvar = 3;
            obj.varFieldIndex = [1 2 11];
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, mesh3d );
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( mesh2d, mesh3d );
            %> time interval
            obj.dt = 0.005;
            %             obj.Cf{1} = 0.0025/1000;
            obj.Cf{1} = 0*ones(size(mesh2d.x));
        end
        
        matTimeSteppingNew( obj );
        
        EntropyAndEnergyCalculation(obj);
        
        AnalysisResult2d( obj );
        AnalysisResult3d( obj );
        
        function NonhydroPostprocess(obj)
            PostProcess = NdgPostProcess(obj.meshUnion(1).mesh2d,strcat('Result/',mfilename,'/2d/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Eta = zeros( Ntime,1 );
            exactEta = zeros( Ntime,1 );
            x0 = 17.5;
            h = obj.H0;
            a = obj.A;
            c = sqrt( obj.gra*obj.Lambda/2/pi*tanh(2*pi*h/obj.Lambda) );
            T = obj.Lambda / c;
            for t = 1:Ntime
                %                 exactEta(t) = -obj.A * cos(2*pi*x0/obj.Lambda + 2*pi*outputTime(t) /obj.T);
                exactEta(t) = obj.A * cos(2*pi/obj.Lambda*x0)*cos(2*pi/T*outputTime(t));
                tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  x0, 0.2, x0, t )-h;
                Eta(t) = tempdata(1);
            end
            figure;
            set(gcf,'position',[50,50,1050,400]);
            plot(outputTime,Eta,'k','LineWidth',1.5);
            hold on;
            set(gca,'YLIM',[-1.3*a, 1.3*a],'Fontsize',15);
            xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
            ylabel({'$\eta\;\rm{(m)}$'},'Interpreter','latex');
            
            %             str = strcat('Hydro',num2str(obj.d),'.fig');
            %             h = openfig(str,'reuse'); % open figure
            %             D1=get(gca,'Children'); %get the handle of the line object
            %             XData1=get(D1,'XData'); %get the x data
            %             YData1=get(D1,'YData'); %get the y data
            %             close(h);
            %             plot(XData1, YData1,'k--','LineWidth',1.5);
            plot(outputTime(1:10:end),exactEta(1:10:end),'ro','markersize',4.5);
            h = legend('SWE','Theory');
            set(h,'fontsize',15);
            legend('boxoff');
        end
        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) =  obj.A * cos(2*pi*obj.mesh2d(m).x/obj.Lambda) - fphys2d{m}(:, :, 4);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            %             option('EddyViscosityType') = enumEddyViscosity.Constant;
            %             option('GOTMSetupFile') = obj.GotmFile;
            %             option('equationType') = enumDiscreteEquation.Strong;
            %             option('integralType') = enumDiscreteIntegral.GaussQuadrature;
            %             option('outputType') = enumOutputFile.VTK;
            %             option('ConstantEddyViscosityValue') = 0;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 1;
            
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedDepth, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, ceil(obj.ChWidth/(obj.ChLength/M)), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

