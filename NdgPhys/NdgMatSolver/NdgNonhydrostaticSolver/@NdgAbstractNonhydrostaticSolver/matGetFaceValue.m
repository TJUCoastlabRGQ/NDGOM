function [fm, fp] = matGetFaceValue( obj, fm, fp, ftype )
%> @brief Function to take the local and adjacent face value
%> @details Function to take the local and adjacent face value
%> @param[in] fm the local face value
%> @param[in] fp the adjacent face value
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface 
%> @param[out] fm the local face value with wet and dry interface considered
%> @param[out] fp the adjacent face value with wet and dry interface considered
[fm, fp] = mxGetFaceValue(fm, fp, obj.NonhydroFmPoint, obj.NonhydroFpPoint, ftype);
end