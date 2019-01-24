function NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
%MATCALCULATENONHYDROSTATICPRESSURE 此处显示有关此函数的摘要
%   此处显示详细说明
 NonhydrostaticPressure = StiffMatrix\ NonhydrostaticRHS(:);
%  NewNonhydrostaticPressure = full(StiffMatrix)\ NonhydrostaticRHS(:);
%  data = NonhydrostaticPressure - NewNonhydrostaticPressure;
%  display(max(max(abs(data))));
%  NonhydrostaticPressure = mxEvaluateNonhydrostaticPressure(StiffMatrix, NonhydrostaticRHS(:) );
end

