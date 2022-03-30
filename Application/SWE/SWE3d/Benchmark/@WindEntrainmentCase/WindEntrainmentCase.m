classdef WindEntrainmentCase < SWEBaroclinic3d
    %WINDENTRAINMENTCASE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ChLength = 2000
        
        ChWidth = 2000
        
        H0 = 50
        
        finalTime = 30*3600
        
        GotmFile = fullfile([pwd,'/Application/SWE/SWE3d/Benchmark/@WindEntrainmentCase/gotmturb.nml'])
    end
    
    properties
        MLDfid
    end
    
    properties ( Constant )
        hcrit = 0.001;
    end
    
    methods
        function obj = WindEntrainmentCase(N, Nz, M, Mz)
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            
            obj.Nfield = 17;
            
            obj.fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy','w', 'hw','hc','rho', 'hT', 'hS', 'Tke', 'Eps'};
            
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            obj.SurfBoundNewmannDate(:,:,1) = 0.1027/obj.rho0 * ones(size(obj.mesh2d.x));%0.1
            
            filename = strcat(['Result\WindEntrainmentCase\MixedLayerDepth-',num2str(Mz),'.txt']);
            
            obj.MLDfid = fopen(filename,'w');
        end
        
        function PostProcess( obj )
            close all;
            figure;
            hold on;
            file = dir('Result\WindEntrainmentCase\*.txt');
            str = cell(1);
            MarkStyle = {'r','g','b','m'};
            for i = 1:numel(file)
                data = importdata(['Result\WindEntrainmentCase\',file(i).name]);
                str{1,i} = strcat([file(i).name(isstrprop(file(i).name,'digit')),' Layers']);
                plot(data(:,1)/3600,data(:,2),MarkStyle{i},'Linewidth',1.5);
            end
            data = importdata(['Result\WindEntrainmentCase\',file(1).name]);
            plot(data(:,1)/3600,1.05*0.01*sqrt(data(:,1)/0.01),'k--','Linewidth',1.5);
            legend({str{1,:},'Theory'},'Location','Northwest');
            box on;
            set(gca,'Linewidth',1.2);
            set(gca,'Fontsize',12);
            xlabel('$t(h)$','Interpreter','Latex','Fontsize',12);
            ylabel('$MDL(m)$','Interpreter','Latex','Fontsize',12);
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
                %> For salinity
                fphys{m}(:,:,15) = obj.H0*obj.S0;
                %> For temperature
                Ttop = 20;
                fphys{m}(:,:,14) = obj.H0 * (Ttop - (0.01)^2*obj.H0*obj.rho0/obj.alphaT/obj.gra*(0-mesh3d.z));
%                 fphys{m}(:,:,14) = obj.H0 * obj.T0;
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0;
            end
        end
        
        function matUpdateOutputResult( obj, time, fphys2d, fphys3d )
            matUpdateOutputResult@SWEAbstract3d( obj, time, fphys2d, fphys3d );
            Index = fphys3d{1}(:,:,16) > 1e-5;
            depth = -1 * obj.H0 * min(min(obj.meshUnion.z(Index)));
            fprintf(obj.MLDfid, '%12.8f  %12.8f\n', time, depth );
        end
        
        function matUpdateFinalResult( obj, time, fphys2d, fphys )
            matUpdateFinalResult@SWEAbstract3d( obj, time, fphys2d, fphys );
            fclose(obj.MLDfid);
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
            %doing nothing
        end
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 800;
            option('startTime') = 0.0;
            option('finalTime') = obj.finalTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = obj.finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('EosType') = enumEOSType.Linear;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.0001;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('PhysicalSurfaceRoughnessLength') = 0.003;
            option('PhysicalBottomRoughnessLength') = 0.003;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.None;
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], 1*[ -obj.ChLength/2, obj.ChLength/2 ], ceil(obj.ChLength/M), 1*ceil(obj.ChLength/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end