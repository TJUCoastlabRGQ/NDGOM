function [fluxNHy, fluxNHvy] = matGetUpwindedNumFluxTermY( obj, Edge, hm, hp, hum, hup, hvm, hvp )
%> @brief Function to calculate the upwinded flux term in direction y
%> @details Function to calculate the upwinded flux term in direction y
%> @param[in] Edge the edge object, which can be boundary edge or inner edge
%> @param[in] hm the local water depth value, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[in] hp the adjacent water depth value, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[in] hvm the local flux rate in y direction, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[in] hvp the adjacent flux rate in y direction, with boundary condition at both the boundary and the wet-dry interface considered
%> @param[out] fluxNHy the calculated numerical flux for water depth, h, in the y direction
%> @param[out] fluxNHvy the calculated numerical flux for momentum, hu, in the x direction

Lnorm = Edge.nx .* hum + Edge.ny .* hvm;     Anorm = -Edge.nx .* hup -Edge.ny .* hvp;

% fluxNHy = Edge.ny .* ( hm + hp )./2;fluxNHvy = Edge.ny .* ( hvm + hvp )./2;
% fluxNHy = zeros(size(Edge.ny));fluxNHvy = zeros(size(Edge.ny));
fluxNHy = Edge.ny .* ( hm + hp )./2;fluxNHvy = Edge.ny .* ( hvm + hvp )./2; %Central for water depth, upwind for velocity term

index = (sign(Lnorm) == 1 & sign(Anorm) ~= 1);
% fluxNHy(index) = Edge.ny(index) .* hm(index);
fluxNHvy(index) = Edge.ny(index) .* hvm(index);
% fluxNHy(index) = Edge.ny(index) .* hm(index);

index = (sign(Lnorm) ~= 1 & sign(Anorm) == 1);
% fluxNHy(index) = -Edge.ny(index) .* hp(index);
fluxNHvy(index) = -Edge.ny(index) .* hvp(index);
% fluxNHy(index) = -Edge.ny(index) .* hp(index);

% index = (sign(Lnorm) .* sign(Anorm) == -1);
% fluxNHy(index) = Edge.ny(index) .* (hm(index).*(sign(Lnorm(index)) + 1)*0.5 + hp(index).*(1 - sign(Lnorm(index)) )*0.5);
% fluxNHvy(index) = Edge.ny(index) .* (hvm(index).*(sign(Lnorm(index)) + 1)*0.5 + hvp(index).*(1 - sign(Lnorm(index)))*0.5);

% index = (sign(Lnorm) .* sign(Anorm) == 1);
% fluxNHy(index) = Edge.ny(index) .* (hm(index) + hp(index))*0.5;
% fluxNHvy(index) = Edge.ny(index) .* (hvm(index) + hvp(index))*0.5;
% 
% index = (sign(Lnorm) == 1 & sign(Anorm) == 0);
% fluxNHy(index) = Edge.ny(index) .* (hm(index).*(sign(Lnorm(index)) + 1)*0.5 );
% fluxNHvy(index) = Edge.ny(index) .* (hvm(index).*(sign(Lnorm(index)) + 1)*0.5);
% 
% index = (sign(Lnorm) == 0 & sign(Anorm) == 1);
% fluxNHy(index) = - Edge.ny(index) .* hp(index).*(1 + sign(Anorm(index)) )*0.5;
% fluxNHvy(index) = - Edge.ny(index) .* hvp(index).*(1 + sign(Anorm(index)) )*0.5;
end