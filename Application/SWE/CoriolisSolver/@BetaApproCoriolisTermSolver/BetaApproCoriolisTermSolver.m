classdef BetaApproCoriolisTermSolver < AbstractCoriolisTermSolver
    properties
        f0
        beta
    end
    
    methods
        function obj = BetaApproCoriolisTermSolver( phys, m, n)
            obj = obj@AbstractCoriolisTermSolver( phys.meshUnion );
            obj.f0 = m;
            obj.beta = n;
        end
        
        function evaluateCoriolisTermRHS( obj, physClass, fphys )

            a = obj.f0;
            b = obj.beta;

            for m = 1:physClass.Nmesh 
                
                mesh = physClass.meshUnion(m);
                ind = (mesh.EToR == int8(enumSWERegion.Wet));

    
                Np = physClass.meshUnion(m).cell.Np;
                k = physClass.meshUnion(m).K;
                s = ones(Np,k);
                q = a*s;%f0
                
                % frhs = frhs + (f+by)hv
                physClass.frhs{m}(:,ind,obj.huIndex) = physClass.frhs{m}(:,ind,obj.huIndex)...
                    + (q(:,ind)+b*mesh.y(:,ind)).*(fphys{m}(:,ind,obj.hvIndex));
                
                % frhs = frhs - (f+by)hu
                physClass.frhs{m}(:,ind,obj.hvIndex) = physClass.frhs{m}(:,ind,obj.hvIndex)...
                    - (q(:,ind)+b*mesh.y(:,ind)).*(fphys{m}(:,ind,obj.huIndex));
                
            end
        end
    end
    
end

