function [ node_val ] = project_vert2node( obj, vert_val )
%PROJECT_VERT2NODE Summary of this function goes here
%   Detailed explanation goes here

t1 = 0.5*( (1-obj.t)*vert_val(1,:) + (1+obj.t)*vert_val(5,:) );
t2 = 0.5*( (1-obj.t)*vert_val(2,:) + (1+obj.t)*vert_val(6,:) );
t3 = 0.5*( (1-obj.t)*vert_val(3,:) + (1+obj.t)*vert_val(7,:) );
t4 = 0.5*( (1-obj.t)*vert_val(4,:) + (1+obj.t)*vert_val(8,:) );

% node_val = 0.5*( -(obj.r+obj.s) .* t1 + (1+obj.r) .* t2 + (1+obj.s) .* t3 );
node_val = bsxfun(@times,(1-obj.r).*(1-obj.s)/4,t1) + bsxfun(@times,(1+obj.r).*(1-obj.s)/4,t2)...
    + bsxfun(@times,(1+obj.r).*(1+obj.s)/4,t3)+ bsxfun(@times,(1-obj.r).*(1+obj.s)/4,t4);
% for i = 1:size(t1,2)
%     node_val(:,i) = 0.5*( -(obj.r+obj.s) .* t1(:,i) + (1+obj.r) .* t2(:,i) + (1+obj.s) .* t3(:,i) );
% end
end