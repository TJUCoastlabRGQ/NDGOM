classdef EllipticProblemAssembleRightHandSideTest < EllipticProblemAssembleFinalStiffMatrix3d
    %ELLIPTICPROBLEMASSEMBLERIGHTHANDSIDETEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        matRHS
    end
    
    methods
        function obj = EllipticProblemAssembleRightHandSideTest(N, Nz, M, Mz)
            obj = obj@EllipticProblemAssembleFinalStiffMatrix3d(N, Nz, M, Mz);
        end
        
        
        RHS = matAssembleRHS( obj, fphys, deltatime, PUPX, PUPS, PVPY, PVPS, PWPS, PSPX, PSPY );
        
        function EllipticProblemSolve(obj)
            obj.fphys{1} = rand(obj.meshUnion.cell.Np, obj.meshUnion.K, obj.Nfield);
            deltatime = rand(1);
            PUPX = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PUPS = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PVPY = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PVPS = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PWPS = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PSPX = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PSPY = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            RHS = obj.NonhydrostaticSolver.TestAssembleRHS( obj, obj.fphys, deltatime, PUPX, PUPS, PVPY, PVPS, PWPS, PSPX, PSPY);
            obj.matRHS = obj.matAssembleRHS( obj.fphys, deltatime, PUPX, PUPS, PVPY, PVPS, PWPS, PSPX, PSPY );
            disp('======================For right hand side===========================')
            disp('The maximum value is:')
            disp(max(max(RHS)))
            disp('The maximum difference is:')
            disp(max(max(RHS - obj.matRHS(:))))
            disp('The minimum difference is:')
            disp(min(min(RHS(:) - obj.matRHS(:))))
            disp('====================End right hand side==============================')

        end
    end
end

