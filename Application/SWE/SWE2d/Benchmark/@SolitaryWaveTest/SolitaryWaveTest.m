classdef SolitaryWaveTest < SWEPreBlanaced2d
    %NONHYDROSTATICSOLITARYWAVERUNUPWALL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 1e-3
%         n = 0.01.^2
        n = 0
        rho = 1000
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
%         x0 = 2.7096
%         x0 = 2.705
        
        
        x0 = 2
        Eta0
        U0
        W0
    end
    
    methods
        function obj = SolitaryWaveTest(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.SolitaryWaveRunUpWall(mesh); 
            obj.initPhysFromOptions( mesh );
            [ obj.Ng, obj.xg, obj.yg ] = setGaugePoint();
        end
        
        function Postprocess(obj)
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            Visual = makeVisualizationFromNdgPhys(obj);
            PostProcess.drawAnimation(Visual, 1, 10, 'SmallSolitaryWaveTest');            
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
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
                
                fphys{m}(:,:,1) = d + obj.Eta0 - fphys{m}(:,:,4) ;
                fphys{m}(:,:,2) = fphys{m}(:,:,1) .* obj.U0;
                fphys{m}(:,:,6) = fphys{m}(:,:,1).*obj.W0./2;
                
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
            ftime = 4;
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
    mesh = makeUniformTriMesh(N, [0, 10], [0, deltax], 10/deltax, deltax/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [0, 10], [0, deltax], 10/deltax, deltax/deltax, bctype);
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