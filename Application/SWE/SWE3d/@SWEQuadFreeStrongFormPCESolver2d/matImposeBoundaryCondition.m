function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, fm, fp, fext2d )
%MATIMPOSEBOUNDARYCONDITION 此处显示有关此函数的摘要
%   此处显示详细说明
ind = ( edge.ftype == enumBoundaryCondition.ClampedDepth );
fp(:, ind, 1) = fext2d(:, ind, 1);

ind = ( edge.ftype == enumBoundaryCondition.ClampedVel );
fp(:, ind, 2) = fext2d(:, ind, 2);
fp(:, ind, 3) = fext2d(:, ind, 3);

ind = ( edge.ftype == enumBoundaryCondition.Clamped );
fp(:, ind, 1) = fext2d(:, ind, 3);
fp(:, ind, 2) = fext2d(:, ind, 1);
fp(:, ind, 3) = fext2d(:, ind, 2);

ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
Hun =  fm( :, ind, 2 ) .* nx(:, ind) + fm( :, ind, 3).* ny(:, ind);
Hvn = -fm( :, ind, 2 ) .* ny(:, ind) + fm( :, ind, 3).* nx(:, ind);

fp(:, ind, 2) = - Hun .* nx(:, ind) - Hvn .* ny(:, ind);
fp(:, ind, 3) = - Hun .* ny(:, ind) + Hvn .* nx(:, ind);

end

