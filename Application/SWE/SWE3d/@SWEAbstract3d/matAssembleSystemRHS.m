function SystemRHS = matAssembleSystemRHS( obj, Tempfphys, SystemRHS, EXa, IMa)
%> @brief Function to calculate the right hand of the system corresponding to the discretization of the vertical diffusion term
%> @details
%> Function to calculate the right hand of the system corresponding to the discretization of the vertical diffusion term
%> @param[in] Tempfphys The physical field at the last time step
%> @param[in] SystemRHS The right hand side to be calculated and returned
%> @param[in] EXa The coefficient corresponding to the explicit part
%> @param[in] IMa The coefficient corresponding to the implicit part
%> @param[out] SystemRHS The right hand side calculated already
EXStage = numel(EXa);
IMStage = numel(IMa);
for i = 1:obj.Nvar
    SystemRHS(:,:,i) = Tempfphys(:,:,i);
for j = 1:numel(IMa)
     SystemRHS(:,:,i) = SystemRHS(:,:,i) + EXa(j)*obj.ExplicitRHS3d(:,:,(i-1)*EXStage+j) + ...
         IMa(j)*obj.ImplicitRHS3d(:,:,(i-1)*IMStage+j);
end
     SystemRHS(:,:,i) = SystemRHS(:,:,i) + EXa(EXStage) * obj.ExplicitRHS3d(:,:,(i-1)*EXStage+EXStage);
end 