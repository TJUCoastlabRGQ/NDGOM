classdef SolitaryWavetransformOverAnSurbmergedBar < SWEPreBlanaced2d
    %NONHYDROSTATICSOLITARYWAVERUNUPWALL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
%         gra = 9.81
        
%         n = 0.01.^2
        n = 0
        rho = 1000
        Depth = 0.0762
    end
    
    properties
        dt
        Ng
        xg
        yg
        A = 0.009144
%         A = 0.0085
%         x0 = 5.95
%         x0 = 2.7096
%         x0 = 2.705
        
        
        x0 = -0.8
        Eta0
        U0
        W0
    end
    
    methods
        function obj = SolitaryWavetransformOverAnSurbmergedBar(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;
            obj.SolitaryWaveTransform(mesh); 
            obj.initPhysFromOptions( mesh );
            [ obj.Ng, obj.xg, obj.yg ] = setGaugePoint();
            obj.matSolve;
            obj.outputFile(1).mergeOutputResult;
            obj.Postprocess;
            %             [ obj.cellId, obj.Vg ] = accessGaugePointMatrix( obj );
        end
        
        function Postprocess(obj)
            index = [1 2 3 4];
            Coor = [1.5936 2.7620 3.6510 4.5146];
            d = [0.0762 0.0381 0.0381 0.0381];
            
            XLIM = {[-1,1],[0.5,2.5],[2,4],[3,5]};
            TitlePosition = {[-0.75,0.2],[0.75,0.2],[2.25,0.2],[3.25,0.2]};
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Eta = zeros(numel(Coor), Ntime);
            for ind = 1:numel(index)
                for t = 1:Ntime
                    tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  Coor(ind), 0.005,  Coor(ind), t );
                    Eta(ind, t) = tempdata(1) - d(ind);
                end
            end
%             fpath = 'D:\PhdResearch\Application\SWE\SWE2d\Benchmark\@NonhydrostaticSolitaryWaveRunUpWall\Data';
           fpath ='D:\PhdResearch\Application\SWE\SWE2d\Benchmark\@SolitaryWavetransformOverAnSurbmergedBar\data\';
            for i = 1:numel(Coor)
                str = strcat(num2str(index(i)),'theory','.csv');
                pathstr = strcat(fpath,str);
                str1 = strcat(num2str(index(i)),'measured','.csv');
                tpathstr = strcat(fpath,str1);
                titlestr = ['Positon',blanks(1),num2str(index(i))];
                data = xlsread(pathstr);
                tdata = xlsread(tpathstr);
                figure;
                 set(gcf,'position',[50,50,1050,400]);
                plot(outputTime - 2.62,Eta(i,:)./obj.Depth,'LineWidth',1.5);
                hold on;
                plot(data(:,1),data(:,2),'rs','markersize',2);
                plot(tdata(:,1),tdata(:,2),'r--','LineWidth',1.5);
                xlim(XLIM{i});
                ylim([-0.1,0.25]);
%                 ylabel({'$\eta$',num2str('/'),'$h_0$'},'Interpreter','latex');
                 ylabel({'$\eta/h_0$'},'Interpreter','latex');   
                xlabel({'$t \;\rm {(s)}$'},'Interpreter','latex'); 
                set(gca,'Fontsize',15);
                title(titlestr,'position',TitlePosition{i},'Fontsize',15);
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
                index = ( mesh.x <= 2);
                bot(index) = -0.0762;
                index = (mesh.x >= 2 & mesh.x <= 2.762);
                bot(index) = -0.0762 + (mesh.x(index) - 2)/20;
                index = (mesh.x >= 2.762);
                bot(index) = -0.0381;
                
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
                
                fphys{m}(:,:,1) = obj.Eta0 - fphys{m}(:,:,4) ;
                fphys{m}(:,:,2) = fphys{m}(:,:,1) .* obj.U0;
                fphys{m}(:,:,6) = fphys{m}(:,:,1) .* obj.W0./2;
                
%                 fphys{m}(:,:,1) = d + a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x - (5.95) - c*0) )).^2 - fphys{m}(:,:,4) ;
%                 fphys{m}(:,:,2) = fphys{m}(:,:,1) .* (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x - (5.95) - c * 0) )).^2)./...
%                     (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x - (5.95) - c * 0) )).^2 +d )*c;
            end
        end
        
        function matUpdateExternalField(obj, tloc, ~)
            a = obj.A;
            d = obj.Depth;
            c = sqrt(obj.gra*(a + d));

            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                edge = obj.meshUnion(m).BoundaryEdge;
                nodeid = bsxfun( @plus, edge.FToN1, (edge.FToE(1, :) - 1) .* mesh.cell.Np);
                obj.fext{m}( :, :, 4 ) = obj.fphys{m}( nodeid + mesh.K * mesh.cell.Np * 3 );
                obj.fext{1}(:,:,1) = d + a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - obj.x0 - c*tloc) )).^2 - obj.fext{m}( :, :, 4 );
                obj.fext{1}(:,:,2) = obj.fext{m}(:,:,1) .* (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - obj.x0 - c * tloc) )).^2)./...
                    (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - obj.x0 - c * tloc) )).^2 +d )*c;
            end
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 9;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
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
    mesh = makeUniformQuadMesh(N, [-2, 6], [0, deltax], 8/deltax, deltax/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

function [Ng, xg, yg] = setGaugePoint( )
% gauge points 5 7 8 9
xg = [1.5936 2.7620 3.6510 4.5146];

yg = [0.0025 0.0025 0.0025 0.0025];

Ng = numel(xg);
end