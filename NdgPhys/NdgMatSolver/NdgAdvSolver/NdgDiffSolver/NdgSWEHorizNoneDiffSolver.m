classdef NdgSWEHorizNoneDiffSolver
    %NDGSWEHORIZNONEDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj =  NdgSWEHorizNoneDiffSolver(~)
            
        end
        
        function  matEvaluateDiffRHS(obj, ~, ~)
            %doing nothing
        end
        
        function matClearGlobalMemory(obj)
            %Doing nothing
        end
    end
    
end

