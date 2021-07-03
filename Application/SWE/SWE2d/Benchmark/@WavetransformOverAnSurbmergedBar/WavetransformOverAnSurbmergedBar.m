classdef WavetransformOverAnSurbmergedBar < SWEPreBlanaced2d
    %WAVETRANSFORMOVERANELLIPTICALSHOAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
%         gra = 9.81
        
        rho = 1000
        amplitude = 0.01
        d = 0.4
        T = 2.02
        ChLength = 35
%         ChWidth = 0.05
        ChWidth = 0.1
    end
    
    properties
        initial_fphys
        length
        k
        spgLength = 10; %> sponge region size
        distance %> distance to boundary nodes
        sigma %> sponge strength
        maxSigma %> maximum sponge strength
        SpongeCoefficient
        Ylim = [0 0.1]
        Xlim = [0 30]
    end
    
    methods (Access = public)
        function obj = WavetransformOverAnSurbmergedBar(N, deltax, cellType)
            
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;
            [ mesh ] = obj.makeUniformMesh(N, deltax, cellType);
            obj.initPhysFromOptions( mesh );
            obj.WaveCharacterEstimate;
            obj.outputFieldOrder2d = [1, 2, 3, 6];
                   
            bp = obj.Xlim(2) - obj.spgLength;
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
        
        function ResetPhys(obj)
            obj.fphys = obj.setInitialField;
        end
        
        function VideoPostprocess(obj)
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat('Result/WavetransformOverAnSurbmergedBar/2d','/','WavetransformOverAnSurbmergedBar'));
            %             Ntime = PostProcess.Nt;
            %             outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Visual = makeVisualizationFromNdgPhys(obj);
            PostProcess.drawAnimation(Visual, 1, 6, 'WavetransformOverAnSurbmergedBar', obj.fphys{1}(:,:,4) );
        end
        
        function VisualPostprocess(obj)
           time = 40;
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
            Eta =  obj.amplitude * sin(2*pi/obj.T*time);
       
            
            % Stelling and Zijlema, 2003
            omega = 2*pi/obj.T;
            obj.fext{1}( :, :, 2 ) =  omega*obj.amplitude/obj.k/obj.d*0.5*(1 + tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (Eta+obj.d);
            
        end
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            matEvaluateTopographySourceTerm@SWEPreBlanaced2d( obj, fphys );

            for m = 1:obj.Nmesh
                obj.frhs{m}(:,:,2) = obj.frhs{m}(:,:,2)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,2);
                obj.frhs{m}(:,:,3) = obj.frhs{m}(:,:,3)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,3);
                obj.frhs{m}(:,:,4) = obj.frhs{m}(:,:,4)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,6);                
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
        
        function evaluateSpongeCoefficient(obj, xb)
            obj.SpongeCoefficient = zeros(size(obj.meshUnion(1).x));
            ratio = ( obj.meshUnion(1).x - xb )/obj.spgLength;
            Index = (ratio>0 & ratio <= 1/2);
            obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(4*ratio(Index)-1)/2)./( 1-(4*ratio(Index)-1).^2) ) +1 );
            Index = (ratio>1/2 & ratio <= 1);
            obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(3-4*ratio(Index))/2)./( 1-(3-4*ratio(Index)).^2) )+1 );
        end
        
        function fphys = setInitialField( obj )
% %                         alpha = 20/360*2*pi;
                        fphys = cell( 1, 1 );
                        mesh = obj.meshUnion(1);
                        fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                        fphys{1}(:,:,1) = 0.4;
%                         fphys{1}(:,:,4) = - 0.04;
                        
                        index =  ( 6 <= mesh.x & mesh.x <= 12);
                        fphys{1}(index) =  0.4 - ( mesh.x(index) - 6 ) ./ 20;
%                         fphys{1}(index + 3 * numel(mesh.x)) = - fphys{1}(index);
                        
                        index = ( 12 <= mesh.x & mesh.x <= 14);
                        fphys{1}(index) =  0.1;
%                         fphys{1}(index + 3 * numel(mesh.x)) = - fphys{1}(index);
                        
                        index = ( 14 <= mesh.x & mesh.x <= 17);
                        fphys{1}(index) =  0.1 + ( mesh.x(index) - 14 ) ./ 10;
%                         fphys{1}(index + 3 * numel(mesh.x)) = - fphys{1}(index);
                        
                        index = ( 18.95 <= mesh.x & mesh.x <= 23.95);
                        fphys{1}(index) =  0.4 - ( mesh.x(index) - 18.95 ) ./ 25;
                        
                        index = ( 23.95 <= mesh.x );
                        fphys{1}(index) = 0.2;
%                         fphys{1}(index + 3 * numel(mesh.x)) = - fphys{1}(index);
            
                        fphys{1}(:,:,4) = -fphys{1}(:,:,1);

%                         fphys{1}(:,:,1) = 0.4;
%                         fphys{1}(:,:,4) = -1 * fphys{1}(:,:,1);
                        obj.initial_fphys = fphys{1};
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 45;
            outputIntervalNum = 2500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
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
        
        function [ mesh ] = makeUniformMesh(obj, N, deltax, type)
            bctype = [...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.ClampedVel, ...
                enumBoundaryCondition.ZeroGrad];
%                 enumBoundaryCondition.ZeroGrad];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, obj.Xlim, obj.Ylim, ceil((obj.Xlim(2) - obj.Xlim(1))/deltax), ceil((obj.Ylim(2) - obj.Ylim(1))/deltax), bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, obj.Xlim, obj.Ylim, ceil((obj.Xlim(2) - obj.Xlim(1))/deltax), ceil((obj.Ylim(2) - obj.Ylim(1))/deltax), bctype);% 20/0.1 22/0.05  %4/0.025, 1/0.0125,
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func
        
    end
    
end



