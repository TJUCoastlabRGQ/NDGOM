classdef NdgNoneVertDiffSolver
    %NDGNONEVERTDIFFSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        
        function obj = NdgNoneVertDiffSolver( physClass )
            %Doing nothing
        end
        
        function fphys = matUpdateImplicitVerticalDiffusion( obj, ~, SystemRHS, ~, ~, ~, ~, ~, ~, ~, ~, ~)
%             ( physClass, SystemRHS, IMa(intRK,intRK), dt, intRK, Stage );
              fphys = SystemRHS;
        end
    end
    
end

