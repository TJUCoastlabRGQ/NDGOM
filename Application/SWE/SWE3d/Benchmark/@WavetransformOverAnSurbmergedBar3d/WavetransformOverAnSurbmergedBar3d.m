classdef WavetransformOverAnSurbmergedBar3d < SWEBarotropic3d
    %WAVETRANSFORMOVERANELLIPTICALSHOAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        rho = 1000
        amplitude = 0.01
        hcrit = 0.01
        d = 0.4
        AD = 0
        T = 2.02
        ChLength = 35
        %         ChWidth = 0.05
        ChWidth = 0.1
    end
    
    properties
        SurfaceBoundaryEdgeType = 'Dirichlet'
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
        Xlim = [0 35]
    end
    
    methods (Access = public)
        function obj = WavetransformOverAnSurbmergedBar3d( N, Nz, M, Mz )
            % setup mesh domain
            [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [1 2 3 11];
            obj.Nfield = 11;
            obj.Nvar = 3;
            obj.varFieldIndex = [1 2 11];
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, mesh3d );
            
            bp = obj.Xlim(2) - obj.spgLength;
            ind = obj.meshUnion.yc > bp; % right part is sponge region
            obj.meshUnion.EToR(ind) = enumSWERegion.Sponge;
            
            obj.WaveCharacterEstimate;
            
            obj.Limiter = NdgVertLimiter3d(obj.meshUnion(1));
            
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
            PostProcess = NdgPostProcess(obj.mesh2d(1),strcat('Result/WavetransformOverAnSurbmergedBar3d/2d','/','WavetransformOverAnSurbmergedBar3d'));
            %             Ntime = PostProcess.Nt;
            %             outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Visual = makeVisualizationFromNdgPhys(obj);
            PostProcess.drawAnimation(Visual, 1, 6, 'WavetransformOverAnSurbmergedBar', obj.fphys2d{1}(:,:,4) );
%             function drawAnimation( obj, Visual, fildId, frameRate, videoName, topography )
        end
        
        function VisualPostprocess(obj)
            time = 10;
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat('Result/WavetransformOverAnSurbmergedBar3d/2d','/','WavetransformOverAnSurbmergedBar3d'));
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
            set(gca,'Ylim',[0.385,0.415]);
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
        
        function matUpdateExternalField( obj, time, ~ , ~)
            Eta =  obj.amplitude * sin(2*pi/obj.T*time);
            % Stelling and Zijlema, 2003
            omega = 2*pi/obj.T;
            
%             h3d = zeros(size(obj.fext3d{1}(:,:,1)));
%             h2d = zeros(size(obj.fext2d{1}(:,:,1)));
%             Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth & (all(obj.meshUnion(1).BoundaryEdge.xb == obj.Xlim(1))' ));
%             h3d(:,Index) = Eta + obj.d;
%             Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth & (all(obj.meshUnion(1).BoundaryEdge.xb == obj.Xlim(2))' ));
%             h3d(:,Index) = obj.d;
%             obj.fext3d{1}(:,:,3) = h3d;
%             
%             
%             Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth & (all(obj.mesh2d.BoundaryEdge.xb  == obj.Xlim(1)))');
%             h2d(:,Index) = Eta + obj.d;
%             Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth & (all(obj.mesh2d.BoundaryEdge.xb  == obj.Xlim(2)))');
%             h2d(:,Index) = obj.d;            
%             obj.fext2d{1}(:,:,3) = h2d;

            hu3d = zeros(size(obj.fext3d{1}(:,:,1)));
            h3d = zeros(size(obj.fext3d{1}(:,:,1)));
            hu2d = zeros(size(obj.fext2d{1}(:,:,1)));
            h2d = zeros(size(obj.fext2d{1}(:,:,1)));
            Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel);
            hu3d(:,Index) = omega*obj.amplitude/obj.k/(obj.d + obj.AD )*0.5*(1+tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (Eta + obj.d + obj.AD);
            obj.fext3d{1}(:,:,1) = hu3d;
            Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
            h3d(:,Index) = Eta + obj.d + obj.AD;
            obj.fext3d{1}(:,:,3) = h3d;
            
            
            Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel );
            hu2d(:,Index) = omega*obj.amplitude/obj.k/(obj.d + obj.AD )*0.5*(1+tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (Eta + obj.d + obj.AD);
            obj.fext2d{1}(:,:,1) = hu2d;
            Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
            h2d(:,Index) = Eta + obj.d + obj.AD;            
            obj.fext2d{1}(:,:,3) = h2d;
        end
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            matEvaluateTopographySourceTerm@SWEBarotropic3d( obj, fphys );
            
            for m = 1:obj.Nmesh
                obj.frhs{m}(:,:,1) = obj.frhs{m}(:,:,1)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,1);
                obj.frhs{m}(:,:,2) = obj.frhs{m}(:,:,2)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,2);
                obj.frhs{m}(:,:,3) = obj.frhs{m}(:,:,3)...
                    - 10 * obj.SpongeCoefficient.* fphys{m}(:,:,11);
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
        
        function [ fphys2d, fphys ] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
%             fphys = cell( obj.Nmesh, 1 );
%             for m = 1 : obj.Nmesh
%                 mesh2d = obj.mesh2d(m);
%                 mesh3d = obj.meshUnion(m);
%                 fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
%                 fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
%                 % bottom elevation
%                 fphys2d{m}(:, :, 4) =  -obj.d * ones(size(mesh2d.x));
%                 %water depth
%                 fphys2d{m}(:,:,1) = -1 .* fphys2d{m}(:, :, 4);
%                 %                  fphys2d{m}(:,:,1) = 2.89677;
%             end
            %                         alpha = 20/360*2*pi;
                                    fphys = cell( 1, 1 );
                                    mesh = obj.meshUnion(1);
                                    mesh2d = obj.mesh2d(1);
                                    fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                                    fphys2d{1} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                                    fphys2d{1}(:,:,1) = 0.4 + obj.AD;
            %                         fphys{1}(:,:,4) = - 0.04;
            
                                    index =  ( 6 <= mesh2d.x & mesh2d.x <= 12);
                                    fphys2d{1}(index) =  0.4 - ( mesh2d.x(index) - 6 ) ./ 20 + obj.AD;
            %                         fphys{1}(index + 3 * numel(mesh.x)) = - fphys{1}(index);
            
                                    index = ( 12 <= mesh2d.x & mesh2d.x <= 14);
                                    fphys2d{1}(index) =  0.1 + obj.AD;
            %                         fphys{1}(index + 3 * numel(mesh.x)) = - fphys{1}(index);
            
                                    index = ( 14 <= mesh2d.x & mesh2d.x <= 17);
                                    fphys2d{1}(index) =  0.1 + ( mesh2d.x(index) - 14 ) ./ 10 + obj.AD;
            %                         fphys{1}(index + 3 * numel(mesh.x)) = - fphys{1}(index);
            
                                    index = ( 18.95 <= mesh2d.x & mesh2d.x <= 23.95);
                                    fphys2d{1}(index) =  0.4 - ( mesh2d.x(index) - 18.95 ) ./ 25 + obj.AD;
            
                                    index = ( 23.95 <= mesh2d.x );
                                    fphys2d{1}(index) = 0.2 + obj.AD;
            
                                    fphys2d{1}(:,:,4) = -fphys2d{1}(:,:,1);
                                    obj.initial_fphys = fphys{1};
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 40;
            outputIntervalNum = 2500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            %             option('EddyViscosityType') = enumEddyViscosity.Constant;
            %             option('GOTMSetupFile') = obj.GotmFile;
            %             option('equationType') = enumDiscreteEquation.Strong;
            %             option('integralType') = enumDiscreteIntegral.GaussQuadrature;
            %             option('outputType') = enumOutputFile.VTK;
            %             option('ConstantEddyViscosityValue') = 0;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 1;
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
        
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedVel, ...
    enumBoundaryCondition.ZeroGrad ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, ceil(obj.ChWidth/(obj.ChLength/M)), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );

end



