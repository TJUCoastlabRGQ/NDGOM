classdef NdgSWEHorizNoneDiffSolver < NdgHorizDiffSolver
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
    end
    
    methods( Access = protected )
        function  matUpdateViscosity(obj)
            %doing nothing
        end
        %> this function is used to update the penalty parameter adopted in
        %> IP form
        function matUpdatePenaltyParameter(obj)
            %doing nothing
        end
    end
    
end

