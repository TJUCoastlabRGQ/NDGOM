function [ FluxM ] = matEvaluateSurfFlux( obj, edge, nx, ny, fm )
%MATEVALUATESURFFLUX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
FluxM = nx .* fm(:,:,2) + ny .* fm(:,:,3);

end

