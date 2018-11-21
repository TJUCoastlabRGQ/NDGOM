classdef NonhydrostaticSolitaryWaveRunUpWall < SWEConventional2d
    %NONHYDROSTATICSOLITARYWAVERUNUPWALL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 1e-3
        n = 0.01.^2
        rho = 1000
    end
    
    properties
        dt
        Ng
        xg
        yg
    end
    
    methods
        function obj = NonhydrostaticSolitaryWaveRunUpWall(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEConventional2d();
            obj.initPhysFromOptions( mesh );
            [ obj.Ng, obj.xg, obj.yg ] = setGaugePoint();
            %             [ obj.cellId, obj.Vg ] = accessGaugePointMatrix( obj );
        end
        
        function Postprocess(obj)
            d = 0.218;
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Eta5 = zeros(Ntime,1);
            for t = 1:Ntime
                tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  obj.xg(1), obj.yg(1),  obj.yg(1), t );
                Eta5(t) = tempdata(1) - d;
            end
            fpath = 'D:\Research\NDG-FEM\SWE\SWE2d\Benchmark\@NonhydrostaticSolitaryWaveRunUpWall\Data';
            pathstr = strcat(fpath,'\point5.csv');
            data = xlsread(pathstr);
            
            figure;
            plot(outputTime,Eta5,'k','LineWidth',1.5);
            hold on;
            plot(data(:,1),data(:,2)/100,'k--','LineWidth',1.5);
            legend('Simulated','gauged');
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
            a = 0.054;
            d = 0.218;
            c = sqrt(obj.gra*(a + d));
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = zeros(size(mesh.x));
                index = (mesh.x >= 15.04 & mesh.x <= 19.40);
                bot(index) = (mesh.x(index) - 15.04)/53;
                index = (mesh.x >= 19.40 & mesh.x <= 22.32);
                bot(index) = 4.36/53 + (mesh.x(index) - 19.4)/150;
                index = (mesh.x >= 22.32 & mesh.x <= 23.2);
                bot(index) = 4.36/53 + 2.92/150 + (mesh.x(index) - 22.32)/13;
                
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
                
                fphys{m}(:,:,1) = d + a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x - (5.95) - c*0) )).^2 - fphys{m}(:,:,4) ;
                fphys{m}(:,:,2) = fphys{m}(:,:,1) .* (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x - (5.95) - c * 0) )).^2)./...
                    (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x - (5.95) - c * 0) )).^2 +d )*c;
            end
        end
        
        function matUpdateExternalField(obj, tloc, ~)
            a = 0.054;
            d = 0.218;
            c = sqrt(obj.gra*(a + d));
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                edge = obj.meshUnion(m).BoundaryEdge;
                nodeid = bsxfun( @plus, edge.FToN1, (edge.FToE(1, :) - 1) .* mesh.cell.Np);
                obj.fext{m}( :, :, 4 ) = obj.fphys{m}( nodeid + mesh.K * mesh.cell.Np * 3 );
                obj.fext{1}(:,:,1) = d + a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - (5.95) - c*tloc) )).^2 - obj.fext{m}( :, :, 4 );
                obj.fext{1}(:,:,2) = obj.fext{m}(:,:,1) .* (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - (5.95) - c * tloc) )).^2)./...
                    (a*(sech(sqrt(3*a/(4*(d)^3)) * (mesh.x(nodeid) - (5.95) - c * tloc) )).^2 +d )*c;
            end
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 30;
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
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.SlipWall];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 23.2], [0, deltax], 23.2/deltax, deltax/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [0, 23.2], [0, deltax], 23.2/deltax, deltax/deltax, bctype);
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