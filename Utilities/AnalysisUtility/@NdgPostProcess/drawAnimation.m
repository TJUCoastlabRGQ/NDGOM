function drawAnimation( obj, Visual, fildId, frameRate, videoName, topography )

% video = VideoWriter( [videoName, '.avi'] );
% video.FrameRate = frameRate;
% open(video);
varId = fildId;
for t = 1:1:obj.Nt
    field = obj.accessOutputResultAtStepNum( t );
    for m = 1:obj.Nmesh
        Visual.drawResult( field{m}(:,:,varId) + topography );
    end
%     view(0,1);
%     zlim([-0.012, 0.012]);
    zlim([-0.08, 0.08]);
%     zlim([0.217, 0.228]);
%     frame = getframe( gcf );
%     writeVideo( video, frame );
    drawnow;
end
% close( video );


