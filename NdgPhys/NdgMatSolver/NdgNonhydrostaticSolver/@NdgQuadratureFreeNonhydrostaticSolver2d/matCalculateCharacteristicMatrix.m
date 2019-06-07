function [ termX, termY] = matCalculateCharacteristicMatrix(obj, mesh,  BoundaryEdge, InnerEdge, VariableX, VariableY, ftype)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges
[fmy, fpy] = InnerEdge.matEvaluateSurfValue( VariableY );  [fmx, fpx] = InnerEdge.matEvaluateSurfValue( VariableX );      
[fmy, fpy] = obj.matGetFaceValue(fmy, fpy, ftype);  [fmx, fpx] = obj.matGetFaceValue(fmx, fpx, ftype);
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx;
fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;

% fluxSx = ( ( fmx + fpx ) ./ 2 - 1/2 .* ( fmx - fpx ) ) .* InnerEdge.nx;
% fluxSy = ( ( fmy + fpy ) ./ 2 - 1/2 .* ( fmy - fpy ) ) .* InnerEdge.ny;

termX = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMX, fluxPX);
termY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMY, fluxPY);

% termY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxPY, fluxSy);
% termX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);



[fmy, fpy] = BoundaryEdge.matEvaluateSurfValue( VariableY );    [fmx, fpx] = BoundaryEdge.matEvaluateSurfValue( VariableX );     
fpy = obj.matImposeNonhydroRelatedBoundaryCondition(fmy, fpy, ftype, obj.EidBoundaryType);
fpx = obj.matImposeNonhydroRelatedBoundaryCondition(fmx, fpx, ftype, obj.EidBoundaryType);
%% test first


%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;
fluxSY = BoundaryEdge.ny.*(fpy + fmy)./2; fluxSX = BoundaryEdge.nx.*(fpx + fmx)./2;

% fluxSX = ( ( fmx + fpx ) ./ 2 - 1/2 .* ( fmx - fpx ) ) .* BoundaryEdge.nx;
% fluxSY = ( ( fmy + fpy ) ./ 2 - 1/2 .* ( fmy - fpy ) ) .* BoundaryEdge.ny;

termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);
termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

[VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral(mesh, cell2mat(VariableX), cell2mat(VariableY));

termY = termY + VolumeIntegralY;
termX = termX + VolumeIntegralX;
end