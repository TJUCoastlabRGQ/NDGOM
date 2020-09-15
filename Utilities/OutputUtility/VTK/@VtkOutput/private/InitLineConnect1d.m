function [ Ncell, EToV ] = InitLineConnect1d(N)
%InitQuadConnect2d - Description
%
% Syntax: [ Ncell, EToV ] = InitQuadConnect2d(input)
%
% Long description
Ncell = N;
EToV = zeros(2, Ncell);

for sk = 1:Ncell
    EToV(:, sk) = [ sk, sk+1 ]';
end
end