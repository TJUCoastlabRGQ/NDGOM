function NonhydrostaticPressure = matCalculateNonhydrostaticPressure(obj, StiffMatrix, NonhydrostaticRHS)
%MATCALCULATENONHYDROSTATICPRESSURE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
 NonhydrostaticPressure = StiffMatrix\ NonhydrostaticRHS(:);
%  NonhydrostaticPressure = mxEvaluateNonhydrostaticPressure(StiffMatrix, NonhydrostaticRHS(:) );
end

