function drawResult( obj )

timeFrac = 1;
Angle = {[122,66]};

WaveTransformPos = NdgPostProcess(obj.meshUnion(1),strcat('WaveTransformOverAnEllipticalShoal','/','WaveTransformOverAnEllipticalShoal'));
[ time ] = ncread( WaveTransformPos.outputFile{1}, 'time' );
Nt = numel( timeFrac );
timeStep = zeros(Nt, 1);
for n = 1:Nt
    timeSpecific = timeFrac(n) * max(time);
    [ ~, timeStep(n) ] = min( abs(time - timeSpecific) );
end
% draw 3d plot
for n = 1:Nt
    figure;
    
    set(gcf,'Position',[309 198 1252 614]);
    ts = timeStep(n);
    field = WaveTransformPos.accessOutputResultAtStepNum( ts );
    
    % draw 3d plot
    %                 draw3dBottom( obj.meshUnion, obj.fphys{1}(:,:,4) );
    %                 colormap copper;
    eta = (field{1}(:,:,1) + obj.fphys{1}(:,:,4))*100;
    draw3dSurface( obj.meshUnion, eta );
    colormap jet;
    
    %                 xlim([5, 24]);
    set(gca,'Ylim',[-10,12]);
    %                 zlim([ - 5,  10]); %case 0.096
    %                 zlim([ - 5,  15]); %case 0.191
    zlim([ - 5,  10]); %case 0.045
    view(Angle{n});
    %                 set( gca, 'CLim', [(0.32 - 0.02 - 0.32)*100, (0.32 + 0.06 - 0.32)*100] );
    xlabel('$x$ (m)', 'FontSize', 14, 'Interpreter', 'Latex');
    ylabel('$y$ (m)', 'FontSize', 14, 'Interpreter', 'Latex');
    zlabel('$\eta$ (cm)', 'FontSize', 14, 'Interpreter', 'Latex');
    box off;
    titletime = timeFrac(n) * max(time); %case 0.181
    title(['$t = $',num2str(roundn(titletime,-2)),' (s)'],'FontSize', 14, 'Interpreter', 'Latex')
    set(gca,'fontsize',14);
    h = colorbar;
    h.Location = 'south';
    h.Position = [0.1 0.04 0.8 0.025];
    h.FontSize = 14;
    h.Limits = [-3,6];
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
            'Vertices', [mesh.x(:), mesh.y(:), (zvar(:))*100], ...
            'Faces', EToV, ...
            'FaceColor', 'w', ...
            'EdgeColor', [0 0 0], ...
            'FaceVertexCData', zvar(:));
        box on;
        grid on;
    end
end