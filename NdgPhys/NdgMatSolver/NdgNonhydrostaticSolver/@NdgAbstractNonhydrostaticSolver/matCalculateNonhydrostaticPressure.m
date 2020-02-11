function NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
%MATCALCULATENONHYDROSTATICPRESSURE 此处显示有关此函数的摘要
%   此处显示详细说明
%% 195.974867s
 NonhydrostaticPressure = StiffMatrix\ NonhydrostaticRHS(:);
 
 %% 633.265648s
%  [L, U] = lu(StiffMatrix);
%  NonhydrostaticPressure = U\(L\NonhydrostaticRHS(:));

%% 207.475206s
% P = symrcm(StiffMatrix);
% StiffMatrix = StiffMatrix(P,P);
% NonhydrostaticRHS = NonhydrostaticRHS(P);
% [L ,U] = lu(StiffMatrix);
% NonhydrostaticPressure = U\(L\NonhydrostaticRHS(:)); 
% NonhydrostaticPressure(P) =  NonhydrostaticPressure;
%%
% C = chol(StiffMatrix);
% NonhydrostaticPressure = C\(C'\NonhydrostaticRHS(:));
end

