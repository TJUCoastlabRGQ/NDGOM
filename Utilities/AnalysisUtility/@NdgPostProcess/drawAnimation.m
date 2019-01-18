function drawAnimation( obj, Visual, fildId, frameRate, videoName )

video = VideoWriter( [videoName, '.avi'] );
video.FrameRate = frameRate;
open(video);
varId = fildId;
for t = 1:5:obj.Nt
    field = obj.accessOutputResultAtStepNum( t );
    for m = 1:obj.Nmesh
        Visual.drawResult( field{m}(:,:,varId) );
    end
    zlim([0.4, 0.5]);
%     zlim([0.217, 0.228]);
    frame = getframe( gcf );
    writeVideo( video, frame );
    drawnow;
end
close( video );


