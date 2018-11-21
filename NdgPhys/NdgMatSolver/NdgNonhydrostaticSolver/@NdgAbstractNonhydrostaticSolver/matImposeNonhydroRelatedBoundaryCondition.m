function fp = matImposeNonhydroRelatedBoundaryCondition(obj, fm, fp, ftype, EidBoundaryType)
%> @brief Function to get the adjacent value at the boundary
%> @details Function to get the adjacent value at the boundary
%> @param[in] fm the local face value
%> @param[in] fp the adjacent face value
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[in] EidBoundaryType the boundary index used to impose the zero nonhydrostatic pressure at the open boundary and the zero-grad condition for nonhydrostatic related auxialary variable at the wall
%> @param[in] ReverseEidBoundaryType the boundary index used to impose the zero-grad nonhydrostatic pressure at the open boundary and the zero condition for nonhydrostatic related auxialary variable at the wall
%> @param[out] fp the adjacent face value at the boundary

fp = mxImposeNonhydroRelatedBoundaryCondition( fm, fp, int8(ftype), EidBoundaryType);

end