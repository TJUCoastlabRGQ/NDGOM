classdef EllipticProblemUpdateConservativeFinalVelocityTest < EllipticProblemMatrixAssembleTest3dNew
    %ELLIPTICPROBLEMUPDATECONSERVATIVEFINALVELOCITYTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        matfphys
    end
    
    methods
        function obj = EllipticProblemUpdateConservativeFinalVelocityTest(N, Nz, M, Mz)
            %ELLIPTICPROBLEMUPDATECONSERVATIVEFINALVELOCITYTEST 构造此类的实例
            %   此处显示详细说明
            obj = obj@EllipticProblemMatrixAssembleTest3dNew(N, Nz, M, Mz);
        end
        
        function EllipticProblemSolve(obj)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            NonhydroPressure = cell(1);
            NonhydroPressure{1} = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            obj.fphys{1} = rand(obj.meshUnion.cell.Np, obj.meshUnion.K, obj.Nfield);
            PSPX = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PSPY = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            obj.matfphys = obj.fphys{1};
            deltatime = rand(1);
            obj.matfphys = obj.matUpdateConservativeFinalVelocity( obj.matfphys, NonhydroPressure, deltatime, PSPX, PSPY );
            obj.fphys = obj.NonhydrostaticSolver.TestUpdateConservativeFinalVelocity( obj, NonhydroPressure{1}, obj.fphys, deltatime, PSPX, PSPY );
            
            disp('=======================For hu========================')
            disp('The maximum value is:')
            disp(max(max(obj.matfphys(:,:,1))))
            disp('The minimum value is')
            disp(min(min(obj.matfphys(:,:,1))))
            disp('The maximum difference is:')
            disp(max(max(obj.matfphys(:,:,1) - obj.fphys{1}(:,:,1))))
            disp('The minimum difference is:')
            disp(min(min(obj.matfphys(:,:,1) - obj.fphys{1}(:,:,1))))
            disp('======================End hu=========================')
            disp('=======================For hv========================')
            disp('The maximum value is:')
            disp(max(max(obj.matfphys(:,:,2))))
            disp('The minimum value is')
            disp(min(min(obj.matfphys(:,:,2))))
            disp('The maximum difference is:')
            disp(max(max(obj.matfphys(:,:,2) - obj.fphys{1}(:,:,2))))
            disp('The minimum difference is:')
            disp(min(min(obj.matfphys(:,:,2) - obj.fphys{1}(:,:,2))))
            disp('======================End hv=========================')
            disp('=======================For hw========================')
            disp('The maximum value is:')
            disp(max(max(obj.matfphys(:,:,11))))
            disp('The minimum value is')
            disp(min(min(obj.matfphys(:,:,11))))
            disp('The maximum difference is:')
            disp(max(max(obj.matfphys(:,:,11) - obj.fphys{1}(:,:,11))))
            disp('The minimum difference is:')
            disp(min(min(obj.matfphys(:,:,11) - obj.fphys{1}(:,:,11))))
            disp('======================End hw=========================')
        end
    end
end

