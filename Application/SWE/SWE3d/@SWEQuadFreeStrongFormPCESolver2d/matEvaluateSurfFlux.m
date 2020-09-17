function [ FluxM ] = matEvaluateSurfFlux( obj, edge, nx, ny, fm )
%MATEVALUATESURFFLUX 此处显示有关此函数的摘要
%   此处显示详细说明
FluxM = nx .* fm(:,:,2) + ny .* fm(:,:,3);

end

