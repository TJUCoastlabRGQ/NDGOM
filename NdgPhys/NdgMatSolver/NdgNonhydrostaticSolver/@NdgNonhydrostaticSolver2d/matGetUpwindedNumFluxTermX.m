function [fluxNHx, fluxNHux] = matGetUpwindedNumFluxTermX( obj, Edge, hm, hp, hum, hup, hvm, hvp )
%> @brief Function to calculate the upwinded flux term in direction x
%> @details Function to calculate the upwinded flux term in direction x
%> @param[in] Edge the edge object, which can be boundary edge or inner edge
%> @param[in] hm the local water depth value, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[in] hp the adjacent water depth value, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[in] hum the local flux rate in x direction, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[in] hup the adjacent flux rate in x direction, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[out] fluxNHx the calculated numerical flux for water depth, h, in the x direction
%> @param[out] fluxNHux the calculated numerical flux for momentum in x direction, hu, in the x direction

Lnorm = Edge.nx .* hum + Edge.ny .* hvm;     Anorm = -Edge.nx .* hup -Edge.ny .* hvp;


% fluxNHx = Edge.nx .* ( hm + hp )./2;  fluxNHux = Edge.nx .* ( hum + hup )./2;
% fluxNHx = zeros(size(Edge.nx));fluxNHux = zeros(size(Edge.nx));
fluxNHx = Edge.nx .* ( hm + hp )./2;  fluxNHux = Edge.nx .* ( hum + hup )./2;  %Central for water depth, upwind for velocity term

index = (sign(Lnorm) == 1 & sign(Anorm) ~= 1);
% fluxNHx(index) = Edge.nx(index) .* hm(index);
fluxNHux(index) = Edge.nx(index) .* hum(index);
% fluxNHx(index) = Edge.nx(index) .* hm(index);

index = (sign(Lnorm) ~= 1 & sign(Anorm) == 1);
% fluxNHx(index) = -Edge.nx(index) .* hp(index);
fluxNHux(index) = -Edge.nx(index) .* hup(index);
% fluxNHx(index) = -Edge.nx(index) .* hp(index);

% index = (sign(Lnorm) .* sign(Anorm) == 1);
% fluxNHx(index) = Edge.nx(index) .* (hm(index) + hp(index))*0.5;
% fluxNHux(index) = Edge.nx(index) .* (hum(index) + hup(index))*0.5;

% index = (sign(Lnorm) == 1 & sign(Anorm) == 0);
% fluxNHx(index) = Edge.nx(index) .* (hm(index).*(sign(Lnorm(index)) + 1)*0.5 );
% fluxNHux(index) = Edge.nx(index) .* (hum(index).*(sign(Lnorm(index)) + 1)*0.5);

% index = (sign(Lnorm) == 0 & sign(Anorm) == 1);
% fluxNHx(index) = - Edge.nx(index) .* hp(index).*(1 + sign(Anorm(index)) )*0.5;
% fluxNHux(index) = - Edge.nx(index) .* hup(index).*(1 + sign(Anorm(index)) )*0.5;
end