classdef WaveTransformOverAnEllipticalShoal3dTest < SWEBarotropic3d
    %WAVETRANSFORMOVERANELLIPTICALSHOAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        rho = 1000
        %         amplitude = 0.0232
        amplitude = 0.0118/2 %0.0232
        d = 0.45
        T = 1
        hcrit = 0.005
        ChLength = 26
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
        spgLength = 4; %> sponge region size
        distance %> distance to boundary nodes
        sigma %> sponge strength
        maxSigma %> maximum sponge strength
        SpongeCoefficient
        Ylim = [-10 -2]
        Xlim = [-10 9.75]
        
        NonhydroIndex
        %         Ylim = [-10 0]
        %         Xlim = [-10 10]
        %         Xlim = [0 0.2]
    end
    
    methods (Access = public)
        function obj = WaveTransformOverAnEllipticalShoal3dTest( N, Nz, Mz )
            %             gmshElementFile = [ fileparts( mfilename('fullpath') ), '/mesh/EllipticShoalElement.msh' ];
            %             gmshBoundaryFile = [ fileparts( mfilename('fullpath') ), '/mesh/EllipticShoalBoundary.msh' ];
            %             [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, Mz, gmshElementFile, gmshBoundaryFile );
            %             SMSFile = [ fileparts( mfilename('fullpath') ), '/mesh/0920.grd' ];
            %             [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, Mz, SMSFile );
            [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, Mz);
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [1 2 3 11];
            obj.Nfield = 11;
            obj.Nvar = 3;
            obj.varFieldIndex = [1 2 11];
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, mesh3d );
            obj.WaveCharacterEstimate;
            
            
            obj.Limiter = NdgVertLimiter3d(obj.meshUnion(1));
            bp = zeros(1,2);
            bp(1) = obj.Ylim(2) - obj.spgLength;
            ind = obj.meshUnion.yc > bp(1); % right part is sponge region
            obj.meshUnion.EToR(ind) = enumSWERegion.Sponge;
            
            %             bp(2) = obj.spgLength + obj.Ylim(1);
            %             ind = obj.meshUnion.yc < bp(2); % left part is sponge region
            %             obj.meshUnion.EToR(ind) = enumSWERegion.Sponge;
            
            
            
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
            PostProcess = NdgPostProcess(obj.mesh2d(1),strcat('Result/WaveTransformOverAnEllipticalShoal3d/2d','/','WaveTransformOverAnEllipticalShoal3d'));
            %             Ntime = PostProcess.Nt;
            %             outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Visual = makeVisualizationFromNdgPhys(obj);
            PostProcess.drawAnimation(Visual, 1, 6, 'WaveTransformOverAnEllipticalShoal3d', obj.fphys2d{1}(:,:,4) );
        end
        
        function VisualPostprocess(obj)
            time = 35;
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat('Result/WaveTransformOverAnEllipticalShoal3d/2d','/','WaveTransformOverAnEllipticalShoal3d'));
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
        
        function matEvaluateWaveMakerCoefficient( obj )
            deltas = 2;
            beltas = 80/deltas/deltas/obj.length;
            I1 = sqrt(pi/beltas)*exp(0);
            omega = 2*pi/obj.T;
            alpha0 = -0.53*(0.5*(-0.53) + 1);
            alpha1 = alpha0 + 1/3;
            D = 2*obj.amplitude*(omega^2 - alpha1 * obj.gra * (obj.k)^4*h^3)/(omega*obj.k*I1*(1-alpha0*(obj.k*h)^2));
        end
        
        function matUpdateExternalField( obj, time, ~, fphys )
            Eta =  obj.amplitude * sin(2*pi/obj.T*time);
            % Stelling and Zijlema, 2003
            omega = 2*pi/obj.T;
            hv3d = zeros(size(obj.fext3d{1}(:,:,1)));
            hw3d = zeros(size(obj.fext3d{1}(:,:,1)));
            
            %             Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel);
            %             ele = obj.meshUnion(1).BoundaryEdge.FToE(1, Index);
            %             % water depth at the boundary
            %             % the sigma coordinate at the boundary
            %             zb = obj.meshUnion(1).z(:,ele);
            %             zb = zb(obj.meshUnion(1).BoundaryEdge.FToN1(:,Index));
            %             % the z level at the boundary
            %             zb = zb.*(Eta + obj.d) + Eta;
            %             TempBoundNonhydroPressure = -obj.NonhydrostaticSolver.rho * obj.gra*obj.amplitude*sin(omega*time)*(1-cosh(obj.k*(zb+obj.d))./cosh(obj.k*obj.d));
            %             TempBoundNonhydroPressure(numel(TempBoundNonhydroPressure(:,1)) - obj.meshUnion.cell.N:end,:) = 0;
            %             obj.NonhydrostaticSolver.BoundNonhydroPressure(:,Index) = TempBoundNonhydroPressure;
            %             hv3d(:,Index) = omega*obj.amplitude/obj.k/(obj.d )*0.5*(1+tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (Eta + obj.d);
            %             obj.fext3d{1}(:,:,2) = hv3d;
            %             Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
            %             h3d(:,Index) = obj.d + Eta;
            %             obj.fext3d{1}(:,:,3) = h3d;
            
            Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel);
            ele = obj.meshUnion(1).BoundaryEdge.FToE(1, Index);
            % water depth at the boundary
            % the sigma coordinate at the boundary
            zb = obj.meshUnion(1).z(:,ele);
            zb = zb(obj.meshUnion(1).BoundaryEdge.FToN1(:,Index));
            % the z level at the boundary
            %             zb = zb.*(Eta + obj.d) + Eta;
            zb = zb.*obj.d;
            TempBoundNonhydroPressure = obj.NonhydrostaticSolver.rho * obj.gra*obj.amplitude*sin(omega*time)*(cosh(obj.k*(zb+obj.d))./cosh(obj.k*obj.d));
            TempBoundNonhydroPressure(numel(TempBoundNonhydroPressure(:,1)) - obj.meshUnion.cell.N:end,:) = 0;
            obj.NonhydrostaticSolver.BoundNonhydroPressure(:,Index) = 0 * TempBoundNonhydroPressure;
            %             hv3d(:,Index) = omega*obj.amplitude*0.5*(1+tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (Eta + obj.d).*(cosh(obj.k*(zb+obj.d))./sinh(obj.k*obj.d));
            %             hv3d(:,Index) = omega*obj.amplitude*0.5*(1+tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (obj.d).*(cosh(obj.k*(zb+obj.d))./sinh(obj.k*obj.d));
            %% The depth-averaged version, value is averaged from bottom to Surface
            % hv3d(:,Index) =  omega*obj.amplitude/obj.k*sinh(obj.k*(obj.d+Eta))/sinh(obj.k*obj.d)*sin(omega*time)*0.5*(1 + tanh((time-3*obj.T)/obj.T));
            % hw3d(:,Index) =  omega*obj.amplitude/obj.k*(cosh(obj.k*(obj.d + Eta))-1)/sinh(obj.k*obj.d)*cos(omega*time)*0.5*(1 + tanh((time-3*obj.T)/obj.T));
            %             %             obj.NonhydrostaticSolver.BoundNonhydroGrad(:,Index) = -1 * obj.NonhydrostaticSolver.rho*obj.gra*obj.amplitude/obj.d*cos(omega*time)*tanh(obj.k*obj.d).*obj.meshUnion(1).BoundaryEdge.ny(:,Index);
            %             obj.NonhydrostaticSolver.BoundNonhydroGrad(:,Index) = -1 *obj.gra*obj.amplitude/obj.d*cos(omega*time)*tanh(obj.k*obj.d).*obj.meshUnion(1).BoundaryEdge.ny(:,Index);
            
            %% The depth-dependent version 1 28.99  43.8
             hv3d(:,Index) = obj.d * omega*obj.amplitude .* (cosh(obj.k*(zb+obj.d))./sinh(obj.k*obj.d)) * sin(omega*time)*0.5*(1 + tanh((time-3*obj.T)/obj.T));
             hw3d(:,Index) = obj.d * omega*obj.amplitude .* (sinh(obj.k*(zb+obj.d))./sinh(obj.k*obj.d)) * cos(omega*time)*0.5*(1 + tanh((time-3*obj.T)/obj.T));
            %% The depth-dependent version 2 23.39 33.75
%             hv3d(:,Index) = (obj.d + Eta) * omega*obj.amplitude .* (cosh(obj.k*(zb+obj.d))./sinh(obj.k*obj.d)) * sin(omega*time)*0.5*(1 + tanh((time-3*obj.T)/obj.T));
%             hw3d(:,Index) = (obj.d + Eta) * omega*obj.amplitude .* (sinh(obj.k*(zb+obj.d))./sinh(obj.k*obj.d)) * cos(omega*time)*0.5*(1 + tanh((time-3*obj.T)/obj.T));
            obj.fext3d{1}(:,:,2) = hv3d;
            %             obj.fext3d{1}(:,:,3) = obj.d + Eta;
            %             obj.fext3d{1}(:,:,3) = obj.d;
            obj.fext3d{1}(:,:,4) = hw3d;
            
            obj.fext2d{1}(:,:,2) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField( hv3d );
            
            
            %             Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel );
            %             hv2d(:,Index) = omega*obj.amplitude/obj.k/(obj.d )*0.5*(1+tanh((time-3*obj.T)/obj.T))*sin(omega*time) * (Eta + obj.d);
            %             hv2d(:,Index) = omega*obj.amplitude/obj.k/(obj.d )*0.5*(1+tanh((time-3*obj.T)/obj.T))*sin(omega*time) * obj.d;
            %% 28.99
            %             hv2d(:,Index) = omega*obj.amplitude/obj.k*0.5*(1 + tanh((time-3*obj.T)/obj.T))*sinh(obj.k*(obj.d + Eta))/sinh(obj.k*obj.d)*sin(omega*time);
            %% The depth-depent version 2 23.39
            %             hv2d(:,Index) = omega*obj.amplitude/obj.k*0.5*(1 + tanh((time-3*obj.T)/obj.T))*sinh(obj.k*(obj.d + Eta))/sinh(obj.k*obj.d)*sin(omega*time);
            
            %             obj.fext2d{1}(:,:,2) = hv2d;
            
                        h3d = zeros(size(obj.fext3d{1}(:,:,1)));
                        Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
                        h3d(:,Index) = obj.d + Eta * 0.5*(1 + tanh((time-10*obj.T)/obj.T)) * 0 ;
                        obj.fext3d{1}(:,:,3) = h3d ;
                        h2d = zeros(size(obj.fext2d{1}(:,:,1)));
                        Index = ( obj.meshUnion(1).mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
                        h2d(:,Index) = obj.d + Eta* 0.5*(1 + tanh((time-10*obj.T)/obj.T)) * 0;
                        obj.fext2d{1}(:,:,3) = h2d ;
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
            
            ys = 5.125;
            %% The following is the momentum source part
            %             soruce = zeros(obj.meshUnion.cell.Np, obj.meshUnion.K);
            %
            %             D = zeros(obj.meshUnion.cell.Np, obj.meshUnion.K);
            %
            %             deltas = 2;
            %
            %             Index = find(all(obj.meshUnion.y>=ys-deltas * obj.length/2/2 & obj.meshUnion.y <= ys + deltas * obj.length/2/2));
            %
            %             deltay = zeros(obj.meshUnion.cell.Np, obj.meshUnion.K);
            %
            %             deltay(:,Index) = obj.meshUnion.y(:,Index) - ys;
            %
            %             beltas = 80/deltas/deltas/obj.length;
            %
            %             I1 = sqrt(pi/beltas)*exp(0);
            %
            %             omega = 2*pi/obj.T;
            %
            %             alpha0 = -0.53*(0.5*(-0.53) + 1);
            %
            %             alpha1 = alpha0 + 1/3;
            %
            %             h = fphys{1}(:,:,4);
            %
            %             D(:,Index) = 2*obj.amplitude*(omega^2 - alpha1 * obj.gra * (obj.k)^4*h(:,Index).^3)./(omega*obj.k*I1*(1-alpha0*(obj.k*h(:,Index)).^2));
            %
            %             soruce(:,Index) = -obj.gra * 2 * beltas .* deltay(:,Index).*exp(-beltas*deltay(:,Index).^2).*D(:,Index)./omega*sin(-1*omega*time);
            %             if time <= 3*obj.T
            %                 soruce(:,Index) = soruce(:,Index) * (1 - exp(-2*time/obj.T));
            %             end
            %
            %             obj.frhs{1}(:,:,2) = obj.frhs{1}(:,:,2) + soruce;
            %% The following is the water depth source part
            %             soruce = zeros(obj.meshUnion.mesh2d.cell.Np, obj.meshUnion.mesh2d.K);
            %
            %             D = zeros(obj.meshUnion.mesh2d.cell.Np, obj.meshUnion.mesh2d.K);
            %
            %             deltas = 2;
            %
            %             Index = find(all(obj.meshUnion.mesh2d.y>=ys-deltas * obj.length/2/2 & obj.meshUnion.mesh2d.y <= ys + deltas * obj.length/2/2));
            %
            %             deltay = zeros(obj.meshUnion.mesh2d.cell.Np, obj.meshUnion.mesh2d.K);
            %
            %             deltay(:,Index) = obj.meshUnion.mesh2d.y(:,Index) - ys;
            %
            %             beltas = 80/deltas/deltas/obj.length;
            %
            %             I1 = sqrt(pi/beltas)*exp(0);
            %
            %             omega = 2*pi/obj.T;
            %
            %             alpha0 = -0.53*(0.5*(-0.53) + 1);
            %
            %             alpha1 = alpha0 + 1/3;
            %
            %             h = fphys2d{1}(:,:,1);
            %
            %             D(:,Index) = 2*obj.amplitude*(omega^2 - alpha1 * obj.gra * (obj.k)^4*h(:,Index).^3)./(omega*obj.k*I1*(1-alpha0*(obj.k*h(:,Index)).^2));
            %
            %             obj.NonhydroIndex = Index;
            %
            %             soruce(:,Index) = D(:,Index) .* sin(-1*omega*time) .*exp(-beltas*deltay(:,Index).^2);
            %             if time <= 3*obj.T
            %                 soruce(:,Index) = soruce(:,Index) * (1 - exp(-2*time/obj.T));
            %             end
            %
            %             obj.frhs2d{1}(:,:,1) = obj.frhs2d{1}(:,:,1) + soruce;
            
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
            ratio = ( obj.meshUnion(1).y - yb(1) )/obj.spgLength;
            Index = (ratio>0 & ratio <= 1/2);
            obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(4*ratio(Index)-1)/2)./( 1-(4*ratio(Index)-1).^2) ) +1 );
            Index = (ratio>1/2 & ratio <= 1);
            obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(3-4*ratio(Index))/2)./( 1-(3-4*ratio(Index)).^2) )+1 );
            
            %             ratio = ( yb(2) - obj.meshUnion(1).y )/obj.spgLength;
            %             Index = (ratio>0 & ratio <= 1/2);
            %             obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(4*ratio(Index)-1)/2)./( 1-(4*ratio(Index)-1).^2) ) +1 );
            %             Index = (ratio>1/2 & ratio <= 1);
            %             obj.SpongeCoefficient(Index) = 1/4*( tanh( sin(pi*(3-4*ratio(Index))/2)./( 1-(3-4*ratio(Index)).^2) )+1 );
            
        end
        
        function [ fphys2d, fphys ] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( 1, 1 );
            mesh = obj.meshUnion(1);
            mesh2d = obj.mesh2d(1);
            fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            fphys2d{1} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
            
            alpha = 20/360*2*pi;
            tempx = mesh2d.x*cos(alpha)+mesh2d.y*sin(alpha);
            tempy = mesh2d.y*cos(alpha) - mesh2d.x*sin(alpha);
            
            index =  (tempy >= -5.84);
            fphys2d{1}(index) =  max(0.1077, obj.d-0.02*(5.84+tempy(index)));
            fphys2d{1}(~index) =  obj.d;
            
            index = (((tempx/4).^2+(tempy/3).^2)<=1);
            fphys2d{1}(index) = fphys2d{1}(index)+0.3-...
                0.5*sqrt(1-(tempx(index)/5).^2-(tempy(index)/3.75).^2);
            
            %             fphys2d{1}(:,:,4) = -fphys2d{1}(:,:,1);
            fphys2d{1}(:,:,1) = obj.d;
            fphys2d{1}(:,:,4) = -obj.d;
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 100;
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
            option('ConstantHorizontalEddyViscosityValue') = 0.1;
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


% function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, Mz, SMSFile )
%
% [ mesh2d ] = makeSMSFileUMeshUnion2d( N, SMSFile );
%
% cell = StdPrismTri( N, Nz );
% zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
% mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
% mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
% mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
% mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
% mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
% mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
%
% Index = ( mesh3d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth);
% mesh3d.BoundaryEdge.ftype(Index) = enumBoundaryCondition.ClampedVel;
%
% Index = ( mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth);
% mesh2d.BoundaryEdge.ftype(Index) = enumBoundaryCondition.ClampedVel;
%
% % [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% % [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
%
% end

% function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, Mz, gmshFileElement, gmshBoundaryFile )
%
% mesh2d = makeGmshFileUMeshUnion2d( N, gmshFileElement, gmshBoundaryFile );
%
% cell = StdPrismTri( N, Nz );
% zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
% mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
% mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
% mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
% mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
% mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
% mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% % [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% % [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
%
% end

% function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, Mz)
%
% bctype = [ ...
%     enumBoundaryCondition.ClampedVel, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall ];
%
% mesh2d = makeUniformQuadMesh( N, ...
%     [ -10, -5 ], [ -10, 6 ], 20, 80, bctype);
% cell = StdPrismQuad( N, Nz );
% zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
% mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
% mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
% mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
% mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
% mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
% mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% % [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
%
% end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, Mz)

bctype = [ ...
    enumBoundaryCondition.ClampedVel, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -10, -9.9 ], [ -10, -2 ], 1, 8/0.1, bctype);
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
