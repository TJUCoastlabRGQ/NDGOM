function matEvaluateRHS( obj, fphys2d, fphys3d )
%MATEVALUATERHS Summary of this function goes here
%   Detailed explanation goes here

obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d, fphys3d, obj.fext2d);

matEvaluateRHS@NdgPhysMat(obj, fphys3d);

end
