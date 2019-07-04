classdef SolitaryWaveRunOnBeach1d < SWEWD1d
    %WAVETRANSFORMOVERANELLIPTICALSHOAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 0.0001
        d = 1
    end
    
    properties
        Xlim = [-100 20]
        H0
        U0
        W0
        P0
    end
    
    methods (Access = public)
        function obj = SolitaryWaveRunOnBeach1d(N, deltax)
            
            obj = obj@SWEWD1d();
            %             obj.hmin = 1e-3;
            [ mesh ] = obj.makeUniformMesh(N, deltax);
            obj.initPhysFromOptions( mesh );
            obj.outputFieldOrder = [1, 2, 3, 6];
        end
        %> Compared numerical water elevation with measured data
        CheckGaugeResult( obj );
        
        function ResetPhys(obj)
            obj.fphys = obj.setInitialField;
        end
        
        function VideoPostprocess(obj)
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
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
        end
        
    end
    
    methods(Access = protected)
        
        function fphys = setInitialField( obj )
            
            fphys = cell( 1, 1 );
            mesh = obj.meshUnion(1);
            fphys{1} = zeros(mesh.cell.Np, mesh.K, obj.Nfield);
            obj.Solitarywave(mesh);
            
            bot = zeros(size(mesh.x));
            index =  ( mesh.x >= -19.85);
            bot(index) =  (mesh.x(index) + 19.85)./19.85;
            
            fphys{1}(:,:,1) = obj.H0 - bot;
            index = fphys{1}(:,:,1)<=0;
            fphys{1}(index) = 0;
            
            temphu = obj.U0 .* fphys{1}(:,:,1);
            index = temphu<=10^(-4);
            temphu(index) = 0;
            fphys{1}(:,:,2) = temphu;
            
            fphys{1}(:,:,3) = bot;
            %             fphys{1}(:,:,5) = obj.W0 .* fphys{1}(:,:,1);
            temphw = obj.W0 .* fphys{1}(:,:,1);
            index = temphw<=10^(-4);
            temphw(index) = 0;
            fphys{1}(:,:,5) = temphw;
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
            %             option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
        end
        
        
        function [ mesh ] = makeUniformMesh( obj, N, deltax)
            M = ceil(( obj.Xlim(2) - obj.Xlim(1) ) / deltax);
            bcType = [...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.ZeroGrad ];
            [ mesh ] = makeUniformMesh1d( N, obj.Xlim, M, bcType );
        end% func
        
    end
    
end



