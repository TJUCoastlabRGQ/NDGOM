function fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys)
mesh = physClass.meshUnion(1);
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
% NonhydrostaticPressure = reshape(NonhydrostaticPressure, mesh.cell.Np, mesh.K);
NonhydrostaticPressure = zeros(mesh.cell.Np, mesh.K);
NonhydroPre = reshape(NonhydroPre, mesh.cell.Np, obj.WetNum);
for i = 1:numel(obj.WetCellIndex)
    NonhydrostaticPressure(:,obj.WetCellIndex(i)) = NonhydroPre(:,i);
end

NonhydroVolumeflux = 1/2 * NonhydrostaticPressure .* fphys{1}(:,:,1);

NqHx = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);
NqHy = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);

% NqHx = matCalculateUpwindedNonhydroRelatedMatrixX(obj, physClass, BoundaryEdge, InnerEdge,  NonhydroVolumeflux, enumNonhydroBoundaryCondition.Zero);
% NqHy = matCalculateUpwindedNonhydroRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge,  NonhydroVolumeflux, enumNonhydroBoundaryCondition.Zero);


fphys{1}(:,:,6) = fphys{1}(:,:,6) + obj.dt/physClass.rho .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) - obj.dt/physClass.rho .* (NqHx +NonhydrostaticPressure.*obj.bx);
fphys{1}(:,:,3) = fphys{1}(:,:,3) - obj.dt/physClass.rho .* (NqHy +NonhydrostaticPressure.*obj.by);

% UpdataWaterDepth(physClass, mesh, fphys, tempFluxhu, tempFluxhv);

end


function termX = matCalculateUpwindedNonhydroRelatedMatrixX(obj, physClass, BoundaryEdge, InnerEdge, Variable, ftype)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges

mesh = physClass.meshUnion(1);
[fm, fp] = InnerEdge.matEvaluateSurfValue( num2cell(Variable,[1 2]) );       
[fm, fp] = obj.matGetFaceValue(fm, fp, ftype);
%< Inner edge contribution
fluxMX = InnerEdge.nx.*fm;
fluxPX = InnerEdge.nx.*fp; 

% fluxSx = zeros(size(InnerEdge.nx));
fluxSx = InnerEdge.nx.*(fp + fm)./2;   %Central
% fluxSx = InnerEdge.nx.*((fp + fm)./2 + InnerEdge.nx .* (fm - fp) + InnerEdge.ny .* (fm - fp));   %Penalty failed

% Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
% fluxSx(Index) = InnerEdge.nx(Index) .* fm(Index);  %upwinded
% Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
% fluxSx(Index) = -InnerEdge.nx(Index) .* fp(Index);  %upwinded

Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
fluxSx(Index) = InnerEdge.nx(Index) .* fp(Index);  %downwinded
Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
fluxSx(Index) = - InnerEdge.nx(Index) .* fm(Index);  %donwinded

% fluxSx = InnerEdge.nx .* (1 + sign(InnerEdge.nx .* fm))./2 .* fp + InnerEdge.nx .*  (1 + sign( -InnerEdge.nx .* fp))./2 .* fm; % downwind
termX = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMX, fluxPX, fluxSx );

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( num2cell(Variable,[1 2]) );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, ftype, obj.EidBoundaryType);
%% test first
%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fm;

% fluxSx = zeros(size(BoundaryEdge.nx));
fluxSx = BoundaryEdge.nx.*(fp + fm)./2;   %Central
% fluxSx = BoundaryEdge.nx.*((fp + fm)./2 + BoundaryEdge.nx .* (fm - fp) + BoundaryEdge.ny .* (fm - fp));   %Penalty failed

% Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
% fluxSx(Index) = BoundaryEdge.nx(Index) .* fm(Index);  %upwinded
% Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
% fluxSx(Index) = -BoundaryEdge.nx(Index) .* fp(Index);  %upwinded

Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
fluxSx(Index) = BoundaryEdge.nx(Index) .* fp(Index);  %downwinded
Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
fluxSx(Index) = - BoundaryEdge.nx(Index) .* fm(Index);  %donwinded

% fluxSX = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.ny .* fm))./2 .* fm+ BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.ny .* fp))./2 .* fp;
% fluxSX = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.ny .* fm))./2 .* fp - BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.ny .* fp))./2 .* fm;


termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSx);

termX = termX + obj.matGetVolumeIntegralX(mesh, Variable);

end

function termY = matCalculateUpwindedNonhydroRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge, Variable, ftype)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges
mesh = physClass.meshUnion(1);
[fm, fp] = InnerEdge.matEvaluateSurfValue( num2cell(Variable,[1 2]) );       
[fm, fp] = obj.matGetFaceValue(fm, fp, ftype);
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fm;
fluxPY = InnerEdge.ny.*fp; 

% fluxSy = zeros(size(InnerEdge.ny)); 
fluxSy = InnerEdge.ny.*(fp + fm)./2;   %Central
% fluxSy = InnerEdge.ny.*((fp + fm)./2 + InnerEdge.nx .* (fm - fp) + InnerEdge.ny .* (fm - fp));   %Penalty failed


% Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
% fluxSy(Index) = InnerEdge.ny(Index) .* fm(Index);  %upwinded
% Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
% fluxSy(Index) = -InnerEdge.ny(Index) .* fp(Index);  %upwinded

Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
fluxSy(Index) = InnerEdge.ny(Index) .* fp(Index);  %downwinded
Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
fluxSy(Index) = - InnerEdge.ny(Index) .* fm(Index);  %downwinded

% fluxSy = InnerEdge.ny .* (1 + sign(InnerEdge.ny .* fm))./2 .* fm - InnerEdge.ny .*  (1 + sign( -InnerEdge.ny .* fp))./2 .* fp;
% fluxSy = InnerEdge.ny .* (1 + sign(InnerEdge.ny .* fm))./2 .* fp + InnerEdge.ny .*  (1 + sign( -InnerEdge.ny .* fp))./2 .* fm;

termY = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMY, fluxPY, fluxSy );

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( num2cell(Variable,[1 2]) );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, ftype, obj.EidBoundaryType);
%% test first
%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fm;

% fluxSy = zeros(size(BoundaryEdge.ny)); 
fluxSy = BoundaryEdge.ny.*(fp + fm)./2;   %Central
% fluxSy = BoundaryEdge.ny.*((fp + fm)./2 + BoundaryEdge.nx .* (fm - fp) + BoundaryEdge.ny .* (fm - fp));   %Penalty failed


% Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
% fluxSy(Index) = BoundaryEdge.ny(Index) .* fm(Index);  %upwinded
% Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
% fluxSy(Index) =  -BoundaryEdge.ny(Index) .* fp(Index);  %upwinded

Index = ( sign(fm) == 1 & sign(fp) ~= 1 );
fluxSy(Index) = BoundaryEdge.ny(Index) .* fp(Index);  %downwinded
Index = ( sign(fm) ~= 1 & sign(fp) == 1 );
fluxSy(Index) = - BoundaryEdge.ny(Index) .* fm(Index);  %downwinded


termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSy);

termY = termY + obj.matGetVolumeIntegralY(mesh, Variable);

% termY = mesh.ry .* (mesh.cell.Dr * Variable)...
%     + mesh.sy .* (mesh.cell.Ds * Variable);   %local

end
