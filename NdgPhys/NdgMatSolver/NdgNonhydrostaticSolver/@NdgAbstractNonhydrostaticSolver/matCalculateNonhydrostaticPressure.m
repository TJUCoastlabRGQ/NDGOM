function NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
%MATCALCULATENONHYDROSTATICPRESSURE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
 NonhydrostaticPressure = StiffMatrix\ NonhydrostaticRHS(:);
%  NewNonhydrostaticPressure = full(StiffMatrix)\ NonhydrostaticRHS(:);
%  data = NonhydrostaticPressure - NewNonhydrostaticPressure;
%  display(max(max(abs(data))));
%  NonhydrostaticPressure = mxEvaluateNonhydrostaticPressure(StiffMatrix, NonhydrostaticRHS(:) );
end

