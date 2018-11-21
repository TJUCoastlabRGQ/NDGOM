function draw( obj, varargin )
%DRAW Summary of this function goes here
%   Detailed explanation goes here

switch nargin
    case 1 % 
        f = obj.f_Q(:,:,1);
    case 2 % ��������Ϊ�ڶ�����������
        f = varargin{1};
end
if ( isempty(obj.draw_h) || ~isvalid(obj.draw_h))
    % ��ͼ��δ���ƻ򴰿ڱ��ر�
    EToV = ones(obj.mesh.K, 1)*obj.mesh.cell.Fmask(:)';
    EToV = EToV + ( obj.mesh.cell.Np*(0:obj.mesh.K-1) )'...
        *ones(1, obj.mesh.cell.Nfptotal);
    obj.draw_h = patch(...
        'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)], ...
        'Faces', EToV, ...
        'FaceColor', 'interp', ...
        'FaceVertexCData', f(:));
else % ��ͼ�����
    set(obj.draw_h, ...
        'Vertices', [obj.mesh.x(:), obj.mesh.y(:), f(:)],...
        'FaceVertexCData', f(:));
end
end

