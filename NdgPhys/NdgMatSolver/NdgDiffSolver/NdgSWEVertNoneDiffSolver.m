classdef NdgSWEVertNoneDiffSolver
    %NDGSWEVERTNONEDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
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

