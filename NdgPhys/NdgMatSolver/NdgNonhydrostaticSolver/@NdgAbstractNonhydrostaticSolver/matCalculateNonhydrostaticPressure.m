function NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
%MATCALCULATENONHYDROSTATICPRESSURE 此处显示有关此函数的摘要
%   此处显示详细说明
 NonhydrostaticPressure = StiffMatrix\ NonhydrostaticRHS(:);
%  NonhydrostaticPressure = mxEvaluateNonhydrostaticPressure(StiffMatrix, NonhydrostaticRHS(:) );
end

