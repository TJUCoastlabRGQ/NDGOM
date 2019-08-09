classdef SolitaryWaveRunUpConicalIsland < SWEPreBlanaced2d
    %SOLITARYWAVERUNUPCONICALISLAND 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Ratio = 0.096
    end
    
    properties
        H0
        U0
        W0
    end
    
    methods
        function obj = SolitaryWaveRunUpConicalIsland(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-4;
            obj.initPhysFromOptions( mesh );
            
            %             [ obj.Ng, obj.xg, obj.yg ] = setGaugePoint();
            %             [ obj.cellId, obj.Vg ] = accessGaugePointMatrix( obj );
        end
    end
    
    methods(Access = protected)
        
        function fphys = setInitialField( obj )
            Xcenter = 12.96;
            Ycenter = 13.80;
            %             mesh = obj.meshUnion(1);
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.Solitarywave(mesh);
                bot = zeros(size(mesh.x));
                radius = sqrt((mesh.x - Xcenter).^2 + (mesh.y - Ycenter).^2);
                
                index = ( radius  <= 1.1 );
                bot(index) = 0.625;
                
                index = ( radius > 1.1 & radius <= 3.6);
                
                bot(index) = 0.625 - ( radius(index) - 1.1)/4;
                
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
                
                index = ( bot >= obj.H0 );
                
                fphys{m}(~index) = obj.H0(~index);
                
                fphys{m}(:,:,2) = fphys{m}(:,:,1) .* obj.U0;
                fphys{m}(:,:,6) = fphys{m}(:,:,1).*obj.W0;
                
            end
        end
        
        function matUpdateExternalField(obj, tloc, ~)
            
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                edge = obj.meshUnion(m).BoundaryEdge;
                nodeid = bsxfun( @plus, edge.FToN1, (edge.FToE(1, :) - 1) .* mesh.cell.Np);
                obj.fext{m}( :, :, 4 ) = obj.fphys{m}( nodeid + mesh.K * mesh.cell.Np * 3 );
                
          
                d = 0.32;
                A = d * obj.Ratio;
                l = d*sqrt((A+d)/A);
                C0 = l/d*sqrt(obj.gra * d^3/(l^2-d^2));
                h = d + A * ( sech( ( 0 -C0*tloc)/l ) )^2;
                U = C0 * (1 - d/h);
                W = -( A * C0 * d )/( l * h ) * sech( (0 - C0 * tloc)/l )*( -sinh((0 - C0*tloc)/l)/(l*cosh((0 - C0*tloc)/l)^2) * l  );
                
                
                obj.fext{1}(:,:,1) = h;
                obj.fext{1}(:,:,2) = h * U;
                obj.fext{1}(:,:,6) = h * W;
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
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
%             option('FrictionType') = enumSWEFriction.Manning;
%             option('FrictionCoefficient_n') = obj.n;
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
    mesh = makeUniformTriMesh(N, [0, 25], [0, 30], 25/deltax, 30/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [0, 25], [13.5, 14.5], 25/deltax, 1/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func

function [Ng, xg, yg] = setGaugePoint( )
% gauge points 3  6 9 16 22
% we note that the x coordinate of point 3 need to be adjusted to match the
% measured data at the point 3, when the data match well, the point will be
% fixed.
xg = [6.71 9.36 10.36 12.96 15.56];
yg = [13.05 13.8 13.8 11.22 13.8];
Ng = numel(xg);
end