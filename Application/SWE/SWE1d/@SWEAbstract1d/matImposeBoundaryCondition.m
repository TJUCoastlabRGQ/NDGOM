function [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, fM, fP, fext )

%Only slip boundary condition considered
for i = 1:numel(edge.ftype)
switch edge.ftype(i)
    case enumBoundaryCondition.ClampedVel
          fP(1,i,1) = fM(1,i,1);
          fP(1,i,2) = fext(1,i,2);
          fP(1,i,5) = fM(1,i,5);
    case enumBoundaryCondition.SlipWall
        fP(1,i,1) = fM(1,i,1);
        fP(1,i,2) = -fM(1,i,2);
        fP(1,i,5) = fM(1,i,5);
    case enumBoundaryCondition.ZeroGrad
        fP(1,i,1) = fM(1,i,1);
        fP(1,i,2) = fM(1,i,2);
        fP(1,i,5) = fM(1,i,5);        
end
end
end