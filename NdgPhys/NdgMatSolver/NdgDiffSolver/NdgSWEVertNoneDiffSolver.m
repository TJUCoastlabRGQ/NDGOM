classdef NdgSWEVertNoneDiffSolver
    %NDGSWEVERTNONEDIFFSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgSWEVertNoneDiffSolver( physClass )
            %Doing nothing
        end
        
        function fphys = matUpdateImplicitVerticalDiffusion( ~, ~, ~, ~, SystemRHS, ~, ~, ~, ~, ~, ~, ~, ~ )
%             ( physClass, SystemRHS, IMa(intRK,intRK), dt, intRK, Stage );
              fphys = SystemRHS;
        end  
        
        function matClearGlobalMemory(obj)
            %doing nothing
        end
    end
    
end

