classdef NonCoriolisTermSolver < AbstractCoriolisTermSolver
    
    methods
        
        function obj = NonCoriolisTermSolver( phys )
            obj = obj@AbstractCoriolisTermSolver( phys.meshUnion );
        end
        
        function evaluateCoriolisTermRHS( obj, physClass, fphys )
            % do nothing ...
        end
    end
    
end

