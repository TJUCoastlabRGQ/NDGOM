function rhs = matAssembleRightHandSide(obj)
x = obj.meshUnion.x;
y = obj.meshUnion.y;
z = obj.meshUnion.z;
rhs = eval(obj.SecondDiffCexact);
% rhs = -2*ones(size(x));
% Temprhs = diag(obj.meshUnion.J(:,1))*obj.meshUnion.cell.M * rhs;
% rhs = Temprhs;
% Tempdata = -pi^2/4*sin(-pi/2*obj.meshUnion.x);
% rhs = zeros(Np*K, 1);
% for i=1:K
%     rhs((i-1)*Np+1:i*Np) = (diag(obj.meshUnion.J(:,i))*obj.meshUnion.cell.M)\sum((diag(obj.meshUnion.J(:,i))*obj.meshUnion.cell.M)*diag(Tempdata(:,i)),2);
% end
end