function drawResult( obj )
          
%             timeFrac = [0.56 0.6 0.645 0.68]; %case 0.096
%             timeFrac = [0.6175 0.66 0.705 0.75]; %case 0.181
            timeFrac = [0.64 0.675 0.74 0.78]; %case 0.045
            Angle = {[-36,55],[-69,38],[41,47],[43,45]};
            
            conicalPos = NdgPostProcess(obj.meshUnion(1),strcat('SolitaryWaveRunUpConicalIsland','/','SolitaryWaveRunUpConicalIsland'));
%             [ time ] = ncread( conicalPos.outputFile{1}, 'time' ) +  21.2 + 2.6;%case 0.181
%             [ time ] = ncread( conicalPos.outputFile{1}, 'time' ) +  22 + 2.6;%case 0.096
            [ time ] = ncread( conicalPos.outputFile{1}, 'time' ) +  25.6;%case 0.045
            Nt = numel( timeFrac );
            timeStep = zeros(Nt, 1);
            for n = 1:Nt
                timeSpecific = timeFrac(n) * max(time);
                [ ~, timeStep(n) ] = min( abs(time - timeSpecific) );
            end
            % draw 3d plot
            for n = 1:Nt
                figure;
                
                set(gcf,'position',[50,50,265,250]);
                ts = timeStep(n);
                field = conicalPos.accessOutputResultAtStepNum( ts );
                
                % draw 3d plot
                draw3dBottom( obj.meshUnion, obj.fphys{1}(:,:,4) );
                eta = (field{1}(:,:,1) + obj.fphys{1}(:,:,4)-0.32)*100;
%                 eta( field{1}(:,:,1) < 1.5e-3 ) = nan; %case 0.096
                eta( field{1}(:,:,1) < 5e-3 ) = nan; %case 0.181
                draw3dSurface( obj.meshUnion, eta );
                colormap jet;
                
                xlim([5, 24]);
                ylim([5, 24]);
%                 zlim([ - 5,  10]); %case 0.096
%                 zlim([ - 5,  15]); %case 0.191
                zlim([ - 5,  5]); %case 0.045
                view(Angle{n});
                set( gca, 'CLim', [(0.32 - 0.02 - 0.32)*100, (0.32 + 0.06 - 0.32)*100] );
                xlabel('$x$ (m)', 'FontSize', 12, 'Interpreter', 'Latex');
                ylabel('$y$ (m)', 'FontSize', 12, 'Interpreter', 'Latex');
                zlabel('$\eta$ (cm)', 'FontSize', 12, 'Interpreter', 'Latex');
                box off;
%                 time = timeFrac(n) * max(time) + 22 + 2.6; %case 0.096
                titletime = timeFrac(n) * max(time); %case 0.181
%                 time = timeFrac(n) * max(time) + 25.6; %case 0.045
                title(['$t = $',num2str(roundn(titletime,-2)),' (s)'],'FontSize', 12, 'Interpreter', 'Latex')
%                 colorbar;
            end
            
                function handle = draw3dSurface( mesh, zvar )
                    EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
                    EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
                    handle = patch(...
                        'Vertices', [mesh.x(:), mesh.y(:), zvar(:)], ...
                        'Faces', EToV, ...
                        'FaceColor', 'interp', ...
                        'EdgeColor', 'none', ...
                        'FaceVertexCData', zvar(:));
                    box on;
                    grid on;
                end

                function handle = draw3dBottom( mesh, zvar )
                    EToV = ones(mesh.K, 1)*mesh.cell.Fmask(:)';
                    EToV = EToV + ( mesh.cell.Np*(0:mesh.K-1) )' * ones(1, mesh.cell.TNfp);
                    handle = patch(...
                        'Vertices', [mesh.x(:), mesh.y(:), (zvar(:)-0.32)*100], ...
                        'Faces', EToV, ...
                        'FaceColor', [.67, .67, .67], ...
                        'EdgeColor', 'none', ...
                        'FaceVertexCData', zvar(:));
                    box on;
                    grid on;
                end
        end