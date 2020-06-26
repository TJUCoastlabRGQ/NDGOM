classdef NdgNoneVertDiffSolver
    %NDGNONEVERTDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
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

