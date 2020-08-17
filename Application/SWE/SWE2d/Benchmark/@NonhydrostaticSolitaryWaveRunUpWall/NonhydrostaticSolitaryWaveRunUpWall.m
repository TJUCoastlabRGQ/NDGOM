classdef NonhydrostaticSolitaryWaveRunUpWall < SWEWDPreBlanaced2d
    %NONHYDROSTATICSOLITARYWAVERUNUPWALL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
%         n = 0.01.^2
        n = 0
        Depth = 0.218
    end
    
    properties
        dt
        Ng
        xg
        yg
        A = 0.054
%         A = 0.0085
%         x0 = 5.95
        x0 = 2.7096
%         x0 = 2.705
        
        
%         x0 = -5
        Eta0
        U0
        W0
    end
    
    methods
        function obj = NonhydrostaticSolitaryWaveRunUpWall(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEWDPreBlanaced2d();
            obj.hmin = 1e-3;
            obj.SolitaryWaveRunUpWall(mesh); 
            obj.initPhysFromOptions( mesh );
            [ obj.Ng, obj.xg, obj.yg ] = setGaugePoint();
%             obj.matSolve;
        end
        
        function Postprocess(obj)
            index = [5 6 7 8 9];
            Coor = [15.04 17.22 19.4 20.86 22.32];
            d = [0.218 0.1769 0.1357 0.1259 0.1162];
%             PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat('NonhydrostaticSolitaryWaveRunUpWall/','2d/','NonhydrostaticSolitaryWaveRunUpWall'));
            [ result ] = PostProcess.interpolateOutputResultToGaugePoint( Coor, 0.01*ones(size(Coor)), Coor );
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );

            fpath ='Application\SWE\SWE2d\Benchmark\@NonhydrostaticSolitaryWaveRunUpWall\Data';
            for i = 1:numel(Coor)
                str = strcat('\Gauge',num2str(index(i)),'.csv');
                pathstr = strcat(fpath,str);
                titlestr = strcat('Gauge',32,num2str(index(i)));
                data = xlsread(pathstr);
                figure;
                set(gcf,'position',[50,50,500,200]);
                plot(outputTime - 6.7,( result(:, 1, i) - d(i) )*100,'k','LineWidth',1.5);
%                 plot(outputTime - 0.1,( result(:, 1, i) - d(i) )*100,'k','LineWidth',1.5);
                hold on;
                plot(data(:,1),data(:,2)*100,'ro','markersize',2);
                ylim([-1.5,10]);
                xlim([0,25]);
%                 ylim([-0.3,2]);
%                 xlim([0,30]);                
                ylabel({'$\eta \;\rm {(cm)}$'},'Interpreter','latex');   
                xlabel({'$t \;\rm {(s)}$'},'Interpreter','latex'); 
                set(gca,'Fontsize',13.5);
                set(gca,'Fontsize',14);
%                 title(titlestr,'position',[4.6,1.4],'Fontsize',14);
                title(titlestr,'position',[4.8,7],'Fontsize',14);
            end
        end
        
    end
    
    methods(Access = protected)
        function [ cellId, Vg ] = accessGaugePointMatrix( obj )
            cellId = cell( obj.Nmesh, 1 );
            Vg = cell( obj.Nmesh, 1 );
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                [ cellId{m}, Vg{m} ] = mesh.accessGaugePointLocation( obj.xg, obj.yg, obj.xg );
            end
        end
        
        function fphys = setInitialField( obj )
            a = obj.A;
            d = obj.Depth;
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = zeros(size(mesh.x));
                index = (mesh.x >= 15.04 & mesh.x <= 19.40);
                bot(index) = (mesh.x(index) - 15.04)/53;
                index = (mesh.x >= 19.40 & mesh.x <= 22.32);
                bot(index) = 4.36/53 + (mesh.x(index) - 19.4)/150;
                index = (mesh.x >= 22.32 & mesh.x <= 23.22);
                bot(index) = 4.36/53 + 2.92/150 + (mesh.x(index) - 22.32)/13;
                
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
                
                fphys{m}(:,:,1) = obj.Eta0 - fphys{m}(:,:,4) ;
                fphys{m}(:,:,2) = fphys{m}(:,:,1) .* obj.U0;
                fphys{m}(:,:,6) = fphys{m}(:,:,1).*obj.W0;
                
            end
        end
        
        function matUpdateExternalField(obj, tloc, ~)
            a = obj.A;
            d = obj.Depth;

            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                edge = obj.meshUnion(m).BoundaryEdge;
                nodeid = bsxfun( @plus, edge.FToN1, (edge.FToE(1, :) - 1) .* mesh.cell.Np);
                obj.fext{m}( :, :, 4 ) = obj.fphys{m}( nodeid + mesh.K * mesh.cell.Np * 3 );
                
%                 obj.fext{1}(:,:,1) = d + a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - obj.x0 - c*tloc) )).^2 - obj.fext{m}( :, :, 4 );
%                 obj.fext{1}(:,:,2) = obj.fext{m}(:,:,1) .* (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - obj.x0 - c * tloc) )).^2)./...
%                     (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - obj.x0 - c * tloc) )).^2 +d )*c;
                x = 0;
                a = obj.A;
                l = d*sqrt((a+d)/a);
                C0 = l/d*sqrt(obj.gra * d^3/(l^2-d^2));
                h = d + a * ( sech( ( x - obj.x0 -C0*tloc)/l ) )^2;
                U = C0 * (1 - d/h);
                W = -( a * C0 * d )/( l * h ) * sech( (x - obj.x0 - C0 * tloc)/l )*( -sinh((x - obj.x0 - C0*tloc)/l)/(l*cosh((x - obj.x0 - C0*tloc)/l)^2) * l  );
                obj.fext{1}(:,:,1) = h;
                obj.fext{1}(:,:,2) = h * U;
                obj.fext{1}(:,:,6) = h * W;                
            end
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 35;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
            option('FrictionType') = enumSWEFriction.Manning;
            option('FrictionCoefficient_n') = obj.n;
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, deltax, type)
bctype = [...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 23.22], [0, deltax], 23.2/deltax, deltax/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [-15, 23.22], [0, deltax], ceil(38.22/deltax), deltax/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

function [Ng, xg, yg] = setGaugePoint( )
% gauge points 5 7 8 9
xg = [15.04 19.4 20.86 22.32];

yg = [0.01 0.01 0.01 0.01];

Ng = numel(xg);
end