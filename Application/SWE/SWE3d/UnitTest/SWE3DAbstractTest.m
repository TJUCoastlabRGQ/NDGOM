classdef SWE3DAbstractTest < handle
    
    properties
        tol = 1e-10
    end
    
    methods
        
        function Assert(obj, Exact, Numerical )
            Ind1 = size(Exact,1);Ind2 = size(Exact,2);Ind3 = size(Exact,3);
            for i = 1:Ind1
                for j = 1:Ind2
                    for k = 1:Ind3
                        assert( abs(Exact(i,j,k)- Numerical(i,j,k)) <= obj.tol );
                    end
                end
            end
            display('Assertation succeed!');
        end
        
    end
end