function rhs = matAssembleRightHandSide( obj )
x = obj.meshUnion.x;
y = obj.meshUnion.y;
z = obj.meshUnion.z;
rhs = eval(obj.SecondDiffCexact);
end