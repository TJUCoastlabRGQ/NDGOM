classdef ManufacturedSolution2d < SWEPreBlanaced2d
    %MANUFACTUREDSOLUTION2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        PostProcess
        timePoint
        ErrNorm2
        Index
        ExactValue
    end
    
    properties
        h
        b
        eta
        u
        v
        ht        
        mhx
        mhy        
        hut
        mhux
        mhuy
        hvt
        mhvx
        mhvy
    end
    
    properties ( Constant )
        %> channel length
        ChLength = 100;
        ChWidth = 100;
        e = 0.001;
        w = 0.01;
%         w = 2*pi/100;
        d = 0.1;
        hcrit = 0.001;
    end
    
    methods
        function obj = ManufacturedSolution2d(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;
            obj.matGetFunction( );
            obj.initPhysFromOptions( mesh );
            %> time interval
            obj.PostProcess = NdgPostProcess(obj.meshUnion(1),...
                strcat('ManufacturedSolution2d/2d','/','ManufacturedSolution2d'));
            obj.ErrNorm2 = cell(3);
            obj.Index = 1; 
            obj.ExactValue = cell(1);
            [obj.ExactValue{1}(:,:,1), obj.ExactValue{1}(:,:,2), obj.ExactValue{1}(:,:,3)] = ...
               obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, obj.ftime); 
        end
        
        function matPlotErrorNorm2(obj)
            close all
            figure(1);
            hold on;
            color = {'r-','k-','b-'};
            LWidth = 1.5;
            FontSize = 15;
            for i = 1:3
                plot(obj.timePoint, obj.ErrNorm2{i},color{i},'LineWidth',LWidth);
            end
            lendstr = {'$h$',...
                '$hu$','$hv$'};
            legend(lendstr,'Interpreter','Latex','FontSize',FontSize);
            xlabel('$time(s)$','Interpreter','Latex');
            ylabel('$L_2 error$','Interpreter','Latex');
            box on;
        end
        
    end
    methods ( Access = protected )
        
        matCalculateExplicitRHSTerm( obj, fphys2d, fphys, Stage, RKIndex, time);
        
        function matGetFunction(obj)
            syms x y z t;
            obj.eta = ( obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) );
            obj.b = - ( 2 - 0.005*( x + y ));
%             obj.b = 0*x + 0*y - 2;
            obj.h = obj.eta - obj.b;
            obj.ht = diff(obj.h,t);
            obj.u = int( obj.h*obj.d.*( z + 1 ).*sin(obj.w.*(x+t)), z, [-1,0] );
            obj.v = int( obj.h*obj.d.*( z + 1 ).*sin(obj.w.*(y+t)), z, [-1,0] );            
%             obj.u = sin(obj.w.*(x+t));
%             obj.v = sin(obj.w.*(y+t));
            obj.mhx = diff( obj.h * obj.u, x);
            obj.mhy = diff( obj.h * obj.v, y);            
            obj.hut = diff( obj.h* obj.u, t);
            obj.mhux = diff( obj.h * obj.u * obj.u + 0.5 * obj.gra * ( obj.h * obj.h - obj.b * obj.b), x);
            obj.mhuy = diff( obj.h * obj.u * obj.v, y);
            obj.hvt = diff( obj.h* obj.v, t);
            obj.mhvx = diff( obj.h * obj.u * obj.v, x);
            obj.mhvy = diff( obj.h * obj.v * obj.v + 0.5 * obj.gra * ( obj.h * obj.h - obj.b * obj.b), y);
        end
        
        
        function matEvaluateError( obj, fphys, time)
            [ H, hu, hv ] = obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, time);
            fext = cell(1);
            fext{1}(:,:,1) = H;
            fext{1}(:,:,2) = hu;
            fext{1}(:,:,3) = hv;
            Tempfphys = cell(1);
            Tempfphys{1}(:,:,1) = fphys(:,:,1);
            Tempfphys{1}(:,:,2) = fphys(:,:,2);
            Tempfphys{1}(:,:,3) = fphys(:,:,3);
            Err2 = obj.PostProcess.evaluateNormErr2( Tempfphys, fext );
            obj.timePoint(obj.Index) = time;
            obj.ErrNorm2{1}(obj.Index)  = Err2(1);
            obj.ErrNorm2{2}(obj.Index)  = Err2(2);
            obj.ErrNorm2{3}(obj.Index)  = Err2(3);
            obj.Index = obj.Index + 1;
        end
        
        %> set initial function
        function [ fphys] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                x = mesh.x; y = mesh.y; t=0;
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                % bottom elevation
                fphys{m}(:, :, 4) = eval(obj.b);
                
                fphys{m}(:,:,5) = ( mesh.rx .* ( mesh.cell.Dr * fphys{m}(:,:,4) ) ) + ...
                    ( mesh.sx .* ( mesh.cell.Ds * fphys{m}(:,:,4) ) );
                
                fphys{m}(:,:,6) = ( mesh.ry .* ( mesh.cell.Dr * fphys{m}(:,:,4) ) ) + ...
                    ( mesh.sy .* ( mesh.cell.Ds * fphys{m}(:,:,4) ) );
                %water depth
                fphys{m}(:,:,1) = eval(obj.h);
                fphys{m}(:,:,2) = eval(obj.h) .* eval(obj.u);        
                fphys{m}(:,:,3) = eval(obj.h) .* eval(obj.v);                
            end
        end
        
        
        function matUpdateExternalField( obj, time, ~, ~ )
            
            [obj.fext{1}( :, :, 1 ), obj.fext{1}( :, :, 2 ), obj.fext{1}( :, :, 3 )] =...
                obj.matGetExactSolution( obj.meshUnion.BoundaryEdge.xb, obj.meshUnion.BoundaryEdge.yb, time);
            
        end
        
        function [ h, hu, hv ] = matGetExactSolution( obj, x, y, t)
            h = eval( obj.h );
            hu = eval( obj.h ) .* eval( obj.u );
            hv = eval( obj.h ) .* eval( obj.v );
        end
        
        
        function [ option ] = setOption( obj, option )
            ftime = 100;
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
    
    methods ( Access = protected )
        
        function  matEvaluateSourceTerm(obj, fphys, time)
            mesh = obj.meshUnion(1);
            matEvaluateSourceTerm@SWEAbstract2d( obj, fphys);
             x = mesh.x; y = mesh.y; t = time;
            obj.frhs{1}(:,:,1) = obj.frhs{1}(:,:,1) + eval( obj.ht ) +...
                eval( obj.mhx ) + eval( obj.mhy );
            
            obj.frhs{1}(:,:,2) = obj.frhs{1}(:,:,2) + eval( obj.hut ) +...
                eval( obj.mhux ) + eval( obj.mhuy )...
                 + obj.gra .* eval( obj.eta ) .* fphys{1}(:,:,5);            
            
            obj.frhs{1}(:,:,3) = obj.frhs{1}(:,:,3) + eval( obj.hvt ) +...
                eval( obj.mhvx ) + eval( obj.mhvy )...
                + obj.gra .* eval( obj.eta ) .* fphys{1}(:,:,6);
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, deltax, type)
bctype = [...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 100], [0, 100], 100/deltax, 100/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N,[0, 100], [0, 100], 100/deltax, 100/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end

% [ mesh ] = ImposePeriodicBoundaryCondition2d(  mesh, 'West-East' );
% [ mesh ] = ImposePeriodicBoundaryCondition2d(  mesh, 'South-North' );

end% func