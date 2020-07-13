classdef ManufacturedSolution3d < SWEBarotropic3d
    %MANUFACTUREDSOLUTION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        PostProcess
        timePoint
        ErrNorm2
        Index
        ExactValue
        lendstr
    end
    
    properties
        h
        b
        ht
        eta
        u
        v
        u2d
        v2d
        Omega
        hut
        mhux
        mhuy
        mhuz
        hvt
        mhvx
        mhvy
        mhvz
        mh2dx
        mh2dy
    end
    
    properties ( Constant )
        %> channel length
        ChLength = 100;
        ChWidth = 100;
        e = 0.001;
%         w = 0.01;
        w = 2*pi/900;
        d = 0.1;
        hcrit = 0.001;
    end
    
    methods
        function obj = ManufacturedSolution3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.matGetFunction();
            obj.PostInit( );
        end
        
        function matPlotErrorNorm2(obj)
            close all
            figure(1);
            hold on;
            color = {'r-','k-','b-','r-'};
            LWidth = 1.5;
            FontSize = 15;
            for i = 1:numel(obj.ErrNorm2)
                plot(obj.timePoint, obj.ErrNorm2{i},color{i},'LineWidth',LWidth);
            end
            legend(obj.lendstr,'Interpreter','Latex','FontSize',FontSize);
            xlabel('$time(s)$','Interpreter','Latex');
            ylabel('$L_2 error$','Interpreter','Latex');
            box on;
        end
        
        function  matPostProcess( obj )
            %This function is used to plot the analytical solution and the
            %numerical solution
            file = dir([mfilename, '\3d\*.csv']);
            figure;
            hold on;
            for i = 1:numel(file)
                date = xlsread([mfilename,'\3d\',file(i).name]);
                [hu,~,vertvelocity,~] = obj.matGetExactSolution( date(:,7), date(:,8), date(:,9), obj.outputFile3d.OutputTime(i));
                plot(date(:,9)',hu);
                plot(date(:,9)',date(:,1)');
            end

        end
    end
    
    methods ( Access = protected )
        
        [ Omega , W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d, time );
        matCalculateExplicitRHSTerm( obj, fphys2d, fphys, Stage, RKIndex, time);
        
        function PostInit(obj)
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3 4];
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            obj.PostProcess = NdgPostProcess(obj.meshUnion(1),...
                strcat('ManufacturedSolution3d/','3d/','ManufacturedSolution3d'));
            obj.ErrNorm2 = cell(4,1);
            obj.Index = 1;
            obj.ExactValue = cell(1);
            [obj.ExactValue{1}(:,:,1), obj.ExactValue{1}(:,:,2), obj.ExactValue{1}(:,:,3), obj.ExactValue{1}(:,:,4)] = ...
                obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, obj.ftime);
            obj.lendstr = {'$hu$',...
                '$hv$','$h$','$\omega$'};
        end
        function matGetFunction(obj)
            syms x y z t;
%             obj.eta = ( obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) );
            obj.eta = ( obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) * 0  );

%             obj.b = - ( 2 - 0.005*( x + y ));
%             obj.b = - ( 2 - 0.005*( x + y )*0);
            obj.b = -( 2 - sin(2*pi/900*x)); 
            obj.h = obj.eta - obj.b;
            obj.u = obj.h*obj.d.*( z + 1 ).*sin(obj.w.*(x+t));
            obj.v = obj.h*obj.d.*( z + 1 ).*sin(obj.w.*(y+t)) * 0;
            obj.ht = diff(obj.h,t);
            obj.u2d = int( obj.u, z, [-1,0] );
            obj.v2d = int( obj.v, z, [-1,0] );
            obj.Omega = int( diff(obj.h * obj.u2d - obj.h * obj.u, x) + diff(obj.h * obj.v2d - obj.h * obj.v, y), z, [-1, z] );
            obj.hut = diff( obj.h* obj.u, t);
            obj.mhux = diff( obj.h * obj.u * obj.u + 0.5 * obj.gra * ( obj.h * obj.h - obj.b * obj.b), x);
            obj.mhuy = diff( obj.h * obj.u * obj.v, y);
            obj.mhuz = diff( obj.u * obj.Omega, z);
            obj.hvt = diff( obj.h* obj.v, t);
            obj.mhvx = diff( obj.h * obj.u * obj.v, x);
            obj.mhvy = diff( obj.h * obj.v * obj.v + 0.5 * obj.gra * ( obj.h * obj.h - obj.b * obj.b), y);
            obj.mhvz = diff( obj.v * obj.Omega, z);
            obj.mh2dx = diff( obj.h * obj.u2d, x);
            obj.mh2dy = diff( obj.h * obj.v2d, y);
            
        end
        
        function matEvaluateError( obj, fphys, time)
            fext = cell(1);
            [ fext{1}(:,:,4), fext{1}(:,:,1), fext{1}(:,:,2), fext{1}(:,:,3)] = ...
                obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time);
            Tempfphys = cell(1);
            Tempfphys{1}(:,:,1) = fphys(:,:,1);
            Tempfphys{1}(:,:,2) = fphys(:,:,2);
            Tempfphys{1}(:,:,3) = fphys(:,:,3);
            Tempfphys{1}(:,:,4) = fphys(:,:,4);
            Err2 = obj.PostProcess.evaluateNormErr2( Tempfphys, fext );
            obj.timePoint(obj.Index) = time;
            obj.ErrNorm2{1}(obj.Index)  = Err2(1);
            obj.ErrNorm2{2}(obj.Index)  = Err2(2);
            obj.ErrNorm2{3}(obj.Index)  = Err2(3);
            obj.ErrNorm2{4}(obj.Index)  = Err2(4);
            obj.Index = obj.Index + 1;
        end
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                % bottom elevation
                x = mesh2d.x; y = mesh2d.y; t = 0;
                
                fphys2d{m}(:, :, 4) = eval(obj.b);
                
                fphys2d{m}(:,:,5) = ( mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                
                fphys2d{m}(:,:,6) = ( mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                %water depth
                fphys2d{m}(:,:,1) = eval( obj.h );
                fphys2d{m}(:, :, 2) = eval( obj.u2d ).* eval( obj.h );
                fphys2d{m}(:, :, 3) = eval( obj.v2d ).* eval( obj.h );
                x = mesh3d.x; y = mesh3d.y; z = mesh3d.z; t = 0;
                fphys{m}(:,:,4) = eval( obj.h );
                fphys{m}(:,:,1) = eval( obj.h ) .* eval( obj.u );
                fphys{m}(:,:,2) = eval( obj.h ) .* eval( obj.v );
                fphys{m}(:,:,3) = eval( obj.Omega );
                %> Z in extended three dimensional fields
                fphys{m}(:,:,6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                %> '$\eta$' in extended three dimensional fields
                fphys{m}(:,:,7) = fphys{m}(:,:,4) + fphys{m}(:,:,6);
                
            end
        end
        
        function  matEvaluateSourceTerm(obj, fphys, time)
            mesh = obj.meshUnion(1);
            mesh2d = obj.mesh2d(1);
            matEvaluateSourceTerm@SWEAbstract3d( obj, fphys);
            x = mesh2d.x; y = mesh2d.y; t = time;
            %[Problematic
            obj.frhs2d{1}(:,:,1) = obj.frhs2d{1}(:,:,1) + eval( obj.ht ) + eval(obj.mh2dx) + eval(obj.mh2dy) ;
            
            x = mesh.x; y = mesh.y; z = mesh.z;
            obj.frhs{1}(:,:,1) = obj.frhs{1}(:,:,1) + eval( obj.hut ) +...
                eval( obj.mhux ) + eval( obj.mhuy )...
                + eval( obj.mhuz ) + obj.gra .* eval( obj.eta ) .* fphys{1}(:,:,8);
            
            
            obj.frhs{1}(:,:,2) = obj.frhs{1}(:,:,2) + eval( obj.hvt ) +...
                eval( obj.mhvx ) + eval( obj.mhvy )...
                + eval( obj.mhvz ) + obj.gra .* eval( obj.eta ) .* fphys{1}(:,:,9);
        end
        
        function matUpdateExternalField( obj, t, ~, ~ )
            x = obj.mesh3d.BoundaryEdge.xb; y = obj.mesh3d.BoundaryEdge.yb;
            z = obj.mesh3d.BoundaryEdge.zb;
            obj.fext3d{1}( :, :, 1 ) = eval( obj.h ) .* eval( obj.u );
            obj.fext3d{1}( :, :, 2 ) = eval( obj.h ) .* eval( obj.v );
            obj.fext3d{1}( :, :, 4 ) = eval( obj.h );
            x = obj.mesh2d.BoundaryEdge.xb; y = obj.mesh2d.BoundaryEdge.yb;
            obj.fext2d{1}( :, :, 1 ) = eval( obj.h );
            x = obj.mesh2d.BoundaryEdge.xb; y = obj.mesh2d.BoundaryEdge.yb;
            obj.fext2d{1}( :, :, 2 ) = eval( obj.h ) .* eval( obj.u2d );
            obj.fext2d{1}( :, :, 3 ) = eval( obj.h ) .* eval( obj.v2d );
            
        end
        
        function  [ hu, hv, Omega, h ] = matGetExactSolution( obj, x, y, z, t)
            hu = eval( obj.h ) .* eval( obj.u );
            hv = eval( obj.h ) .* eval( obj.v );
            h = eval( obj.h );
            Omega = eval( obj.Omega );
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 86.4;
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
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped ];

% mesh2d = makeUniformQuadMesh( N, ...
%     [ 0, obj.ChLength ], [ 0, obj.ChWidth ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);
mesh2d = makeUniformQuadMesh(N,[0, 900], [0, M], 900/M, 1, bctype);
cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zeros(mesh2d.Nv, 1) - ones(mesh2d.Nv, 1);
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end