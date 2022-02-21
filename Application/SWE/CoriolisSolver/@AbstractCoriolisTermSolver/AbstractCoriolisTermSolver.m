classdef AbstractCoriolisTermSolver < handle
    
    properties(Constant)
    end
    
    properties
        huIndex
        
        hvIndex
    end
    
    methods
        function obj = AbstractCoriolisTermSolver( mesh )
            if mesh.type == enumMeshDim.Two
                obj.huIndex = 2;
                obj.hvIndex = 3;
            elseif mesh.type == enumMeshDim.Three
                obj.huIndex = 1;
                obj.hvIndex = 2;
            end
        end
    end
    
    methods(Abstract)
        frhs = evaluateCoriolisTermRHS( physClass, fphys )
    end
    
end

