function matCalculatePenaltyParameter( obj, mesh)
%> @brief Function to calculate the penalty parameter C11 for the LDG flux calculation
%> @details
%> Function to  calculate the penalty parameter C11 for the LDG flux calculation.
%> This function is programmed according to Hesthaven and Warburton 2007.
%> @param[in] mesh The mesh object
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;

end