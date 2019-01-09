classdef WaveTransformOverAnEllipticalShoal < SWEPreBlanaced2d
    %WAVETRANSFORMOVERANELLIPTICALSHOAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 1e-3
        rho = 1000
        amplitude = 0.0232
        d = 0.575
        T = 1
    end
    
    properties
        initial_fphys
        length
        k
    end
    
    methods (Access = public)
        function obj = WaveTransformOverAnEllipticalShoal(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.initPhysFromOptions( mesh );
            obj.WaveCharacterEstimate;
        end
        %> Compared numerical water elevation with measured data
        CheckGaugeResult( obj );
        
        function ResetPhys(obj)
            obj.fphys = obj.setInitialField;
        end

        function Postprocess(obj)  
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Visual = makeVisualizationFromNdgPhys(obj);
            PostProcess.drawAnimation(Visual, 1, 10, 'RK33NHydrostaticWave');
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
        
        function WaveCharacterEstimate(obj)
           f = @(L) L - obj.gra*(obj.T)^2/(2*pi)*tanh(2*pi/L*obj.d);
           obj.length = fzero(f,[1 2]);
           obj.k = 2*pi/obj.length;
        end
        
        function matUpdateExternalField( obj, time, ~ )
            Eta =  obj.amplitude * cos(2*pi*time - pi/2);
%             Eta =  obj.amplitude * cos(-2*pi*time);
            %NHWAVE  surface elevation and velocity
            
%             syms z;
%             obj.fext{1}( :, :, 1 ) = obj.d + Eta;
%             aveV = int(2*Eta*pi/1*cosh(obj.k*(z+obj.d))/sinh(obj.k*obj.d),-obj.d,Eta)/(obj.d + obj.amplitude * cos(-2*pi*time));
%             obj.fext{1}( :, :, 3 ) = aveV * obj.fext{1}( :, :, 1 );


            % SWASH velocity
            obj.fext{1}( :, :, 1 ) = obj.d + Eta;
            obj.fext{1}( :, :, 3 ) =  (obj.d + Eta) .* sqrt(obj.gra ./ (obj.d + Eta)) * Eta;
            % Geoclaw water depth and velocity 
%             obj.fext{1}( :, :, 1 ) = d + obj.amplitude * cos(2*pi*time+pi/2);
%             obj.fext{1}( :, :, 3 ) =  2 * obj.fext{1}( :, :, 1 ) .*  (sqrt(obj.gra * obj.fext{1}( :, :, 1 )) - sqrt(obj.gra * d));
            
            % pressure setting Cui Haiyang
%             obj.fext{1}( :, :, 6 ) = -2 * obj.rho * obj.gra * A * cos(2*pi*time+pi/2)*...
%                 (1-sinh(k*obj.fext{1}( :, :, 1 ))./k./obj.fext{1}( :, :, 1 )./cosh(-1*k*0.45));            
%             index = 1;
%             ftype = repmat(obj.meshUnion(1).BoundaryEdge.ftype,[1 obj.meshUnion(1).BoundaryEdge.Nfp])';
%             for i = 1:obj.meshUnion(1).BoundaryEdge.Ne
%                 for j  = 1:obj.meshUnion(1).BoundaryEdge.Nfp
%                     if ftype(index) ~= enumBoundaryCondition.ClampedDepth
%                         obj.fext{1}( j, i, 6 ) = 0;
%                     end
%                     index = index+1;
%                 end
%             end
        end
        
        function fphys = setInitialField( obj )
%             alpha = 0/360*2*pi;
%             fphys = cell( 1, 1 );
%             mesh = obj.meshUnion(1);
%             fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
%             tempx = mesh.x*cos(alpha)+mesh.y*sin(alpha);
%             tempy = mesh.y*cos(alpha) - mesh.x*sin(alpha);
%             index =  (tempy >= -5.84);
%             fphys{1}(index) =  obj.d-0.02*(5.84+tempy(index));
%             fphys{1}(~index) =  obj.d;
%             
%             index = (((tempx/4).^2+(tempy/3).^2)<1);
%             fphys{1}(index) = fphys{1}(index)+0.3-...
%                 0.5*sqrt(1-(tempx(index)/5).^2-(tempy(index)/3.75).^2);
%             fphys{1}(:,:,4) = -fphys{1}(:,:,1);
%  
%             
%             obj.initial_fphys = fphys{1};
              
              fphys = cell( 1, 1 );
              mesh = obj.meshUnion(1);
              fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
              fphys{1}(:,:,1) = obj.d;
              fphys{1}(:,:,4) = -fphys{1}(:,:,1);
            
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 4.5;
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
            option('NumFluxType') = enumSWENumFlux.LF;
%             enumSWENumFlux.LF
%             option('nonhydrostaticType') = enumNonhydrostaticType.Hydrostatic;
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, ~, type)
bctype = [...
    enumBoundaryCondition.ClampedDepth, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [-10, 10], [-10, 12], 20/0.1, 22/0.05, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [-9.2, -9.18], [0, 10], 0.02/0.02, 10/0.02, bctype);% 20/0.1 22/0.05  %4/0.025, 1/0.0125,
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func