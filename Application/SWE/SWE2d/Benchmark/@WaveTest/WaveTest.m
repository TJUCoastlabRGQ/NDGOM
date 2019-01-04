classdef WaveTest < SWEConventional2d
    %WAVETEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 1e-3
        rho = 1000
    end
    
    properties
        dt
    end
    
    methods
        function obj = WaveTest(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEConventional2d();
            obj.initPhysFromOptions( mesh );
        end
        
        function Postprocess(obj)  
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Visual = makeVisualizationFromNdgPhys(obj);
            PostProcess.drawAnimation(Visual, 1, 10, 'NonhydrostaticWaveTest');
            Eta = zeros( Ntime,1 );
            exactEta = zeros( Ntime,1 );
            Lambda = 20;
            x0 = 17.5;
            h = 7.5;
            a = 0.1;
            c = sqrt(obj.gra*Lambda/2/pi*tanh(2*pi*h/Lambda));
            T = Lambda/c;
            for t = 1:Ntime
                exactEta(t) = a * cos(2*pi/Lambda*x0)*cos(2*pi/T*outputTime(t)) + h;
                tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  x0, 0.2, x0, t );
                Eta(t) = tempdata(1);
            end
            figure;
            plot(outputTime,Eta,'k','LineWidth',1.5);
            hold on;
            plot(outputTime,exactEta,'k--','LineWidth',1.5);
            legend('Simulated','Exact');
            xlabel('time(s)');
            ylabel('Water level(m)');
        end
    end
    
    methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                Lambda = 20;
                bot = -7.5;
                k = 2*pi/Lambda;
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
                fphys{m}(:,:,1) = -fphys{m}(:,:,4);
                
                Eta = zeros(size(fphys{m}(:,:,1)));
                PEta = zeros(size(fphys{m}(:,:,1)));
                
                colIndex = all(mesh.x>=5) & all(mesh.x<=25);
                Eta(:,colIndex) = 0.1 * cos(2*pi*mesh.x(:,colIndex)/Lambda);
                PEta(:,colIndex) = 0.1 * sin(2*pi*mesh.x(:,colIndex)/Lambda);
                
                fphys{m}(:,:,1) =  Eta - fphys{m}(:,:,4);
                
                syms z;
                aveU = zeros(size(mesh.x));
                aveW = zeros(size(mesh.x));
                for i = 1:size(mesh.x,1)
                    for j = 1:size(mesh.x,2)
                      aveU(i,j) = int(2*pi/3.59*Eta(i,j)*cosh(k*(z - bot))/sinh(-k*bot),bot,Eta(i,j))/(-bot + Eta(i,j));
%                       aveW(i,j) = int(2*pi/3.59*PEta(i,j)*sinh(k*(z - bot))/sinh(-k*bot),bot,Eta(i,j))/(-bot + Eta(i,j));
                    end 
                end
                fphys{m}(:,:,2) = aveU .* fphys{m}(:,:,1) + 1 * fphys{m}(:,:,1);
                fphys{m}(:,:,6) = aveW .* fphys{m}(:,:,1);
%                 fphys{m}(:,colIndex,2) =  fphys{m}(:,colIndex,1) * 20/3.59;
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 15;
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
%             option('nonhydrostaticType') = NdgNonhydrostaticType.Hydrostatic;
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, deltax, type)
bctype = [...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 20], [0, deltax], 20/deltax, deltax/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [-30, 100], [0, deltax], 130/deltax, deltax/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
