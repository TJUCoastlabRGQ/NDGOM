function rhs2d = matEvaluate2dBoundaryStressRHS(obj, mesh3d, fphys, rhs2d)
%> @brief Function used to calculate the contribution to the external mode due to the bottom stress and surface wind stress
%>
%> More detailed description.
%>
%> @param mesh3d The three dimensional mesh object
%> @param fphys The three dimensional physical field
%> @param rhs2d The two dimensional right hand side calculated with the two dimensional solver
%>
%> @retval rhs2d The two dimensional right hand side with the bottom stress and surface wind stress considered

for n = 1:obj.Nmesh
    %> surface wind stress term
    rhs2d{n}(:,:,2) = rhs2d{n}(:,:,2)  + obj.Taux{n};
    rhs2d{n}(:,:,3) = rhs2d{n}(:,:,3)  + obj.Tauy{n};
    
    edge = mesh3d(n).BottomBoundaryEdge;
    [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
    u = fm(:,:,1)./fm(:,:,6); v = fm(:,:,2)./fm(:,:,6);
    Velocity = sqrt( u.^2 + v.^2 );
    rhs2d{n}(:,:,2) = rhs2d{n}(:,:,2)  - obj.Cf{n} .* u .* Velocity;
    rhs2d{n}(:,:,3) = rhs2d{n}(:,:,3)  - obj.Cf{n} .* v .* Velocity;    
end
end