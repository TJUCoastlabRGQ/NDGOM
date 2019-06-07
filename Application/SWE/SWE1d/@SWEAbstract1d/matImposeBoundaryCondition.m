function [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, fM, fP, fext )

%Only slip boundary condition considered
fP(:,:,2) = - fM(:,:,2);
end