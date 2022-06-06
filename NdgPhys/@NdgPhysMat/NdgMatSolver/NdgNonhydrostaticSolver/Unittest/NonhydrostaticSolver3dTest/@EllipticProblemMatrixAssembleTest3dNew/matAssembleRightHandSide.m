function rhs = matAssembleRightHandSide( obj )
x = obj.meshUnion.x;
y = obj.meshUnion.y;
z = obj.meshUnion.z;
rhs = eval(obj.SecondDiffCexact);
for i = 1:obj.meshUnion.K
    rhs(:,i) = diag(obj.meshUnion.J(:,i)) * obj.meshUnion.cell.M * rhs(:, i);
end
end