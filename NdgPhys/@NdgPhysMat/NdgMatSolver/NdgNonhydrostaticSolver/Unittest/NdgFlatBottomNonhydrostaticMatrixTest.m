classdef NdgFlatBottomNonhydrostaticMatrixTest < NdgFlatBottomWetDryTestWithWallBoundary
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgFlatBottomNonhydrostaticMatrixTest(N, cellType)
            obj = obj@NdgFlatBottomWetDryTestWithWallBoundary(N, cellType);
        end
        
        function matrixX = getExactMatrixX(obj)
            matrixX = ones(size(obj.meshUnion(1).x));
        end
        
        function matrixY = getExactMatrixY(obj)
            matrixY = zeros(size(obj.meshUnion(1).x));
        end
        
        function matrixX = getExactConservativeVariableRelatedMatrixX(obj)
            matrixX = ones(size(obj.meshUnion(1).x));
        end
        
        function matrixY = getExactConservativeVariableRelatedMatrixY(obj)
            matrixY = zeros(size(obj.meshUnion(1).x));
        end
        
        function ExactStiffMatrix = getGlobalStiffMatrix(obj)
            height = obj.fphys{1}(:,:,1);
            mesh = obj.meshUnion(1);
            ExactStiffMatrix = zeros(numel(mesh.x));
            Np = mesh.cell.Np;K = mesh.K;
            for i = 1:K*Np
                gmat = zeros( Np, K );
                gmat(i) = 1.0/2.0;
                ExactStiffMatrix(i,i) = 2*1*1/obj.rho;
                [ Exactqx, Exactqy ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(height(i) .* gmat,[1 2]), num2cell(height(i) .* gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
                [ Exactq2x, Exactq2y ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(Exactqx,[1 2]), num2cell(Exactqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
                ExactStiffMatrix(:,i) = ExactStiffMatrix(:,i) + 1/obj.rho .* Exactqx(:) - 1/obj.rho .* height(:).*Exactq2x(:)-...
                    1/obj.rho .* height(:).*Exactq2y(:);
            end
        end
        
    end
    
    methods(Access = protected)
        
        function fphys = setInitialField( obj )
            fphys = cell(obj.Nmesh, 1);
            fphys{1}(:,:,4) = zeros(size(obj.meshUnion(1).x));
            fphys{1}(:,:,1) = 1 + obj.meshUnion(1).x;
        end
        
        function [ mesh ] = makeUniformMesh(obj, N, type)
            bctype = [...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, [0, 3], [0, 3], 3, 3, bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, [0, 3], [0, 3], 3, 3, bctype);
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func
    end
    
end

