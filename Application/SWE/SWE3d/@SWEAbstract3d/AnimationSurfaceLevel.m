function AnimationSurfaceLevel( obj )
%ANIMATIONSURFACELEVEL Summary of this function goes here
%   Detailed explanation goes here

visual = Visual2d( obj.mesh2d );
OutstepNum = obj.outputFile.outputStep;
% OutstepNum = 2299;
video = VideoWriter( [obj.outputFile.casename, '/', ...
    obj.outputFile.casename, '.avi'] );
video.FrameRate = 20;
open( video );

% initialize axis with the first step output
figure('Color', 'w');
[ fphys2d, ~ ] = obj.outputFile.readOutputResult( obj, 1);
visual.drawResult( fphys2d(:, :, 1) );
zlim([7.1, 8.2]);
zlabel('$\xi$ (m)', 'Interpreter', 'Latex', 'FontSize', 14);
xlabel('$x$ (m)', 'Interpreter', 'Latex', 'FontSize', 14);
% ylabel('$y$ (m)', 'Interpreter', 'Latex', 'FontSize', 14);

for n = 1 : 2 : OutstepNum
    [ fphys2d, ~ ] = obj.outputFile.readOutputResult( obj, n);
    H = fphys2d(:, :, 1);
    
    visual.drawResult( H );
    frame = getframe( gcf );
    writeVideo( video, frame );
end

close( video );

end

