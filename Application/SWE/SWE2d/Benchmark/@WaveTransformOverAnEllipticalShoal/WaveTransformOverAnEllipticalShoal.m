classdef WaveTransformOverAnEllipticalShoal < SWEPreBlanaced2d
    %WAVETRANSFORMOVERANELLIPTICALSHOAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
%         gra = 9.81
        
        rho = 1000
        amplitude = 0.0232
        d = 0.45
        T = 1
        ChLength = 26
%         ChWidth = 0.05
        ChWidth = 0.1
    end
    
    properties
        initial_fphys
        length
        k
        spgLength = 4; %> sponge region size
        distance %> distance to boundary nodes
        sigma %> sponge strength
        maxSigma %> maximum sponge strength
        SpongeCoefficient
        Ylim = [-10 16]
%         Xlim = [-10 10]
        Xlim = [0 0.2]
    end
    
    methods (Access = public)
        function obj = WaveTransformOverAnEllipticalShoal(N, deltax, cellType)
            
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;
            [ mesh ] = obj.makeUniformMesh(N, deltax, cellType);
            obj.initPhysFromOptions( mesh );
            obj.WaveCharacterEstimate;
            obj.outputFieldOrder2d = [1, 2, 3, 6];
                   
            bp = obj.Ylim(2) - obj.spgLength;
            ind = obj.meshUnion.yc > bp; % right part is sponge region
            obj.meshUnion.EToR(ind) = enumSWERegion.Sponge;
            
            %methods from LongXiang Li              
%             Nb = 10;
%             xb = bp * ones( Nb, 1 );
%             yb = linspace( 0, obj.ChWidth, Nb )';
%             obj.evaluateSpongeDistance( xb, yb );
%             dt = matUpdateTimeInterval( obj, obj.fphys );
%             obj.evaluateSpongeStrength( obj.spgLength, 0.05 /dt );
           obj.evaluateSpongeCoefficient(bp);
%            obj.matSolve;
            
        end
        %> Compared numerical water elevation with measured data
        CheckGaugeResult( obj );
        
        drawResult( obj );
        
        function ResetPhys(obj)
            obj.fphys = obj.setInitialField;
        end
        
        function VideoPostprocess(obj)
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
            %             Ntime = PostProcess.Nt;
            %             outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Visual = makeVisualizationFromNdgPhys(obj);
            PostProcess.drawAnimation(Visual, 1, 6, 'WaveTransformOverAnEllipticalShoal', obj.fphys{1}(:,:,4) );
        end
        
        function VisualPostprocess(obj)
           time = 35;
           PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
           outputTime = ncread( PostProcess.outputFile{1}, 'time' );
           [~,Index] = sort(abs(outputTime-time));
           [ fphys ] = PostProcess.accessOutputResultAtStepNum(  Index(1) );
           Visual = makeVisualizationFromNdgPhys( obj );
           Visual.drawResult( fphys{1}(:,:,1)+ obj.fphys{1}(:, :, 4) );
           shading interp; 
%            axis off;
%            hold on;
           xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex','Fontsize',12);
           ylabel({'$y\;\rm{(m)}$'},'Interpreter','latex','Fontsize',12);
           zlabel({'$\eta\;\rm{(m)}$'},'Interpreter','latex','Fontsize',12);
           set(gca,'Ylim',[-10,12]);
           set(gca,'FontSize',12);
           set(gcf,'Position',[309 198 1252 614]);
           view(122, 66);
           h = colorbar;
           h.Location = 'south';
           h.Position = [0.1 0.04 0.8 0.025];
           h.FontSize = 12;
%            box off;
%            grid off;
           
%            Visual.drawResult( obj.fphys{1}(:, :, 4)+fphys{1}(:,:,1) );
%            Visual.drawHandle.FaceColor = [0.5 0.5 0.5];
%            Visual.drawResult( fphys{1}(:,:,1)+ obj.fphys{1}(:, :, 4));
%            Visual.drawResult( obj.fphys{1}(:, :, 4));
%            Visual.drawHandle.FaceColor = [0.5 0.5 0.5];
%            Visual.drawHandle.FaceAlpha = 1;
%            light('Position',[-1 0 0],'Style','local');
        end
        
    end
    
    methods(Access = protected)
        
        function WaveCharacterEstimate(obj)
            f = @(L) L - obj.gra*(obj.T)^2/(2*pi)*tanh(2*pi/L*obj.d);
            obj.length = fzero(f,[1 20]);
            obj.k = 2*pi/obj.length;
        end
        
        function matUpdateExternalField( obj, time, ~ )
%             Eta =  obj.amplitude * cos(2*pi/obj.T*time - pi/2);
            Eta =  obj.amplitude * sin(2*pi/obj.T*time);
            %NHWAVE  surface elevation and velocity
            
            %             syms z;
            %             obj.fext{1}( :, :, 1 ) = obj.d + Eta;
            %             aveV = int(2*Eta*pi/1*cosh(obj.k*(z+obj.d))/sinh(obj.k*obj.d),-obj.d,Eta)/(obj.d + obj.amplitude * cos(-2*pi*time));
            %             obj.fext{1}( :, :, 3 ) = aveV * obj.fext{1}( :, :, 1 );
            % SWASH velocity
%             obj.fext{1}( :, :, 1 ) = obj.d + Eta;
%             obj.fext{1}( :, :, 3 ) =  (obj.d + Eta) .* sqrt(obj.gra ./ (obj.d + Eta)) * Eta;
            %             obj.fext{1}( :, :, 3 ) =  (obj.d + Eta) .* sqrt(obj.gra ./ (obj.d)) * Eta;
            % Geoclaw water depth and velocity
            %             obj.fext{1}( :, :, 1 ) = obj.d + Eta;
            %             obj.fext{1}( :, :, 3 ) =  2 * obj.fext{1}( :, :,
            %             1 ) .*  (sqrt(obj.gra * obj.fext{1}( :, :, 1 )) - sqrt(obj.gra * obj.d));
            
            % Stelling and Zijlema, 2003
            omega = 2*pi/obj.T;
            obj.fext{1}( :, :, 3 ) =  omega*obj.amplitude/obj.k/obj.d*0.5*(1 + tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (Eta+obj.d);
            
%              pressure setting Cui Haiyang
%             obj.fext{1}( :, :, 6 ) = -2 * obj.rho * obj.gra * Eta*...
%                 (1-sinh(obj.k*obj.fext{1}( :, :, 1 ))./obj.k./obj.fext{1}( :, :, 1 )./cosh(obj.k*obj.d)).*obj.fext{1}( :, :, 1 )./2;     %(hq0)/2
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
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            matEvaluateTopographySourceTerm@SWEPreBlanaced2d( obj, fphys );
            
%             for m = 1:obj.Nmesh
%                 obj.frhs{m}(:,:,1) = obj.frhs{m}(:,:,1)...
%                     - obj.sigma.*( fphys{m}(:,:,1) - obj.d );
%                 obj.frhs{m}(:,:,2) = obj.frhs{m}(:,:,2)...
%                     - obj.sigma.*( fphys{m}(:,:,2) - 0 );
%                 obj.frhs{m}(:,:,3) = obj.frhs{m}(:,:,3)...
%                     - obj.sigma.*( fphys{m}(:,:,3) - 0 );
%             end

            for m = 1:obj.Nmesh
                obj.frhs{m}(:,:,2) = obj.frhs{m}(:,:,2)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,2);
                obj.frhs{m}(:,:,3) = obj.frhs{m}(:,:,3)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,3);
            end
        end
        
        %> \brief calculate distance from the boundary
        function evaluateSpongeDistance(obj, xb, yb)
            
            mesh = obj.meshUnion;
            obj.distance = zeros( mesh.cell.Np, mesh.K);
            for i = 1:mesh.K
                if (mesh.EToR(i) ~= enumSWERegion.Sponge)
                    continue;
                end
                
                for n = 1:mesh.cell.Np
                    xi = mesh.x(n, i);
                    yi = mesh.y(n, i);
                    obj.distance(n, i) = min( sqrt( (xi - xb).^2 + (yi - yb).^2 ) );
                end
            end
        end% func        
        
        function evaluateSpongeCoefficient(obj, yb)
            obj.SpongeCoefficient = zeros(size(obj.meshUnion(1).x));
            ratio = ( obj.meshUnion(1).y - yb )/obj.spgLength;
            Index = (ratio>0 & ratio <= 1/2);
            obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(4*ratio(Index)-1)/2)./( 1-(4*ratio(Index)-1).^2) ) +1 );
            Index = (ratio>1/2 & ratio <= 1);
            obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(3-4*ratio(Index))/2)./( 1-(3-4*ratio(Index)).^2) )+1 );
        end
        
        function fphys = setInitialField( obj )
                        alpha = 20/360*2*pi;
                        fphys = cell( 1, 1 );
                        mesh = obj.meshUnion(1);
                        fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                        tempx = mesh.x*cos(alpha)+mesh.y*sin(alpha);
                        tempy = mesh.y*cos(alpha) - mesh.x*sin(alpha);
                        
                        index =  (tempy >= -5.84);
                        fphys{1}(index) =  max(0.1077, obj.d-0.02*(5.84+tempy(index)));
                        fphys{1}(~index) =  obj.d;

%                         index =  (tempy >= -5.484);
%                         fphys{1}(index) =  max(0.1,obj.d-0.02*(5.84+tempy(index)));
%                         fphys{1}(~index) =  obj.d;
                        
%                         index = (mesh.y > 12);
%                         fphys{1}(index) = max(0.1077,obj.d-0.02*(5.84+tempy(index)));
                        
                        index = (((tempx/4).^2+(tempy/3).^2)<=1);
                        fphys{1}(index) = fphys{1}(index)+0.3-...
                            0.5*sqrt(1-(tempx(index)/5).^2-(tempy(index)/3.75).^2);

%                         index = (((tempx/4).^2+(tempy/3).^2)<=1);
%                         fphys{1}(index) = -0.3+...
%                             0.5*sqrt(1-(tempx(index)/5).^2-(tempy(index)/3.75).^2);
                        
                        fphys{1}(:,:,4) = -fphys{1}(:,:,1);
            
                        obj.initial_fphys = fphys{1};
            
%             fphys = cell( 1, 1 );
%             mesh = obj.meshUnion(1);
%             fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
%             fphys{1}(:,:,1) = obj.d;
%             fphys{1}(:,:,4) = -fphys{1}(:,:,1);
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 35;
            outputIntervalNum = 2500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK33;
            option('outputNcfileNum') = 500;            
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
        end
        
        function evaluateSpongeStrength(obj, spongeLength, maxSigma)
            mesh = obj.meshUnion;
            obj.sigma = zeros(mesh.cell.Np, mesh.K);
            p = 3;
            for i = 1:mesh.K
                if mesh.EToR(i) ~= enumSWERegion.Sponge
                    continue;
                end
                obj.sigma(:, i) = maxSigma*abs( obj.distance(:, i)/spongeLength ).^p;
            end
        end% func
        
        function [ mesh ] = makeUniformMesh(obj, N, ~, type)
            bctype = [...
                enumBoundaryCondition.ClampedVel, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, [-10, 10], [-10, 12], 20/0.1, 22/0.05, bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, obj.Xlim, obj.Ylim, ceil((obj.Xlim(2) - obj.Xlim(1))/obj.ChWidth/2), ceil((obj.Ylim(2) - obj.Ylim(1))/obj.ChWidth), bctype);% 20/0.1 22/0.05  %4/0.025, 1/0.0125,
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func
        
    end
    
end

