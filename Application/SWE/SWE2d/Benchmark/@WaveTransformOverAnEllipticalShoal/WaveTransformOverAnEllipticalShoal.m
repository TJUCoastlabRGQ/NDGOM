classdef WaveTransformOverAnEllipticalShoal < SWEConventional2d
    %WAVETRANSFORMOVERANELLIPTICALSHOAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 1e-3
        rho = 1000
        amplitude = 0.0232
    end
    
    properties
        initial_fphys
    end
    
    methods (Access = public)
        function obj = WaveTransformOverAnEllipticalShoal(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEConventional2d();
            obj.initPhysFromOptions( mesh );
        end
        %> Compared numerical water elevation with measured data
        CheckGaugeResult( obj );
    end
    
    methods(Access = protected)
        
        function matUpdateExternalField( obj, time, ~ )
            a = 0.0232;
            obj.fext{1}( :, :, 1 ) = 0.45 + obj.amplitude * cos(2*pi*time);
        end
        
        function fphys = setInitialField( obj )
            alpha = 20/360*2*pi;
            fphys = cell( 1, 1 );
            mesh = obj.meshUnion(1);
            fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            tempx = mesh.x*cos(alpha)+mesh.y*sin(alpha);
            tempy = mesh.y*cos(alpha) - mesh.x*sin(alpha);
            index =  (tempy >= -5.84);
            fphys{1}(index) =  0.45-0.02*(5.84+tempy(index));
            fphys{1}(~index) =  0.45;
            
            index = (((tempx/4).^2+(tempy/3).^2)<1);
            fphys{1}(index) =  0.45-0.02*(5.84+tempy(index))+0.3-...
                0.5*sqrt(1-(tempx(index)/5).^2-(tempy(index)/3.75).^2);
            fphys{1}(:,:,4) = -fphys{1}(:,:,1);
            
            obj.initial_fphys = fphys{1};
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
%             option('nonhydrostaticType') = enumNonhydrostaticType.Hydrostatic;
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, deltax, type)
bctype = [...
    enumBoundaryCondition.ClampedDepth, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [-10, 10], [-10, 12], 20/deltax, 22/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [-10, 10], [-10, 12], 20/deltax, 22/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func