function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, nz, fm, fp, fext3d )

ind = ( edge.ftype == enumBoundaryCondition.Clamped );
fp(:, ind, 1) = fext3d(:, ind, 1);
fp(:, ind, 2) = fext3d(:, ind, 2);
fp(:, ind, 4) = fext3d(:, ind, 4);

ind = ( edge.ftype == enumBoundaryCondition.ClampedDepth );
fp(:, ind, 4) = fext3d(:, ind, 4);

ind = ( edge.ftype == enumBoundaryCondition.ClampedVel );
fp(:, ind, 1) = fext3d(:, ind, 1);
fp(:, ind, 2) = fext3d(:, ind, 2);

ind = ( edge.ftype == enumBoundaryCondition.ZeroGrad );
fp(:, ind, 1) = fm(:, ind, 1);
fp(:, ind, 2) = fm(:, ind, 2);

ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
Hun =  fm( :, ind, 1 ) .* nx(:, ind) + fm( :, ind, 2).* ny(:, ind);
Hvn = -fm( :, ind, 1 ) .* ny(:, ind) + fm( :, ind, 2).* nx(:, ind);

fp(:, ind, 1) = - Hun .* nx(:, ind) - Hvn .* ny(:, ind);
fp(:, ind, 2) = - Hun .* ny(:, ind) + Hvn .* nx(:, ind);
end