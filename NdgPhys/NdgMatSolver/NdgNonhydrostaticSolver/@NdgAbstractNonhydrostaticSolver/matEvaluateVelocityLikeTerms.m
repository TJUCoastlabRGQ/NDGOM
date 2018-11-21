function harmonicityVelocity = matEvaluateVelocityLikeTerms(obj, mesh, fm, fp)
tol = 10^(-8);
index = (fm + fp >= tol);
harmonicityVelocity = zeros(size(fm));
harmonicityVelocity(index) = 2* fm(index) .* fp(index)./(fm(index) + fp(index));
end