classdef NdgFlatBottomWetDryTestWithWallBoundary < NdgNonhydrostaticAbstractTest
    %<@brief Class used to test the nonhydrostatic model considering wetting and drying
    %< In this test, all the boundary are set to be slipwall
    
    methods
        function obj = NdgFlatBottomWetDryTestWithWallBoundary(N, cellType)
            obj = obj@NdgNonhydrostaticAbstractTest(N, cellType);
        end
        
        function EcactWetDryFaceOrder = getExactWetDryFaceOrder(obj)
            %> added on 2019/6/3
            EcactWetDryFaceOrder = [4 5 7 9];
        end
        
        function ExactWetCellIndex=getExactWetCellIndex(obj)
            ExactWetCellIndex = [1 2 3 4 6 7 8 9];
        end
        
        function [ExactZeroFluxBoundary, ExactZeroFluxBoundaryIndex] = getExactZeroFluxBoundary(obj)
            ExactZeroFluxBoundary = [2 3;4 2;6 4;8 1;];
            ExactZeroFluxBoundaryIndex = 4;
        end
        
        function ExactAdjacentDryCellAndFace = getExactAdjacentDryCellAndFace(obj)
            ExactAdjacentDryCellAndFace = [5 1;5 4;5 2;5 3];
        end
        
        function ExactEidBoundaryType = getExactEidBoundaryType(obj)
            ExactEidBoundaryType =  ones(size(obj.meshUnion(1).BoundaryEdge.FToN1));
        end
        
        function ExactReverseEidBoundaryType = getExactReverseEidBoundaryType(obj)
            Nfp = obj.meshUnion(1).cell.Nfp;
            ExactReverseEidBoundaryType = ones(size(obj.meshUnion(1).eidM));
            ExactReverseEidBoundaryType(1:Nfp,1) = -1;
            ExactReverseEidBoundaryType(3*Nfp + 1 :4*Nfp,1) = -1;
            ExactReverseEidBoundaryType(1:Nfp,2) = -1;
            ExactReverseEidBoundaryType(1:2*Nfp,3) = -1;
            ExactReverseEidBoundaryType(3*Nfp + 1 :4*Nfp,4) = -1;
            ExactReverseEidBoundaryType(Nfp + 1 :2*Nfp,6) = -1;
            ExactReverseEidBoundaryType(2*Nfp + 1 :4*Nfp,7) = -1;
            ExactReverseEidBoundaryType(2*Nfp + 1 :3*Nfp,8) = -1;
            ExactReverseEidBoundaryType(Nfp + 1 :3*Nfp,9) = -1;
        end
        
        function ExactGetFaceOuterValue = getExactGetFaceOuterValue(obj)
            Nfp = obj.meshUnion(1).cell.Nfp;
            Nonhydro = obj.fphys{1}(:,:,1);
            ExactGetFaceOuterValue = Nonhydro(obj.meshUnion(1).eidP);
            ExactGetFaceOuterValue(2*Nfp+1:3*Nfp,2) = -1;
            ExactGetFaceOuterValue(Nfp+1:2*Nfp,4) = -1;
            ExactGetFaceOuterValue(3*Nfp+1:4*Nfp,6) = -1;
            ExactGetFaceOuterValue(1:Nfp,8) = -1;
            ExactGetFaceOuterValue(:,5) = 0;
        end
        
        function ExactEidBoundaryType = getUpdatedExactEidBoundaryType(obj)
            Nfp = obj.meshUnion(1).cell.Nfp;
            ExactEidBoundaryType = ones(size(obj.meshUnion(1).eidM));
            ExactEidBoundaryType(2*Nfp+1:3*Nfp,2) = -1;
            ExactEidBoundaryType(Nfp+1:2*Nfp,4) = -1;
            ExactEidBoundaryType(3*Nfp+1:4*Nfp,6) = -1;
            ExactEidBoundaryType(1:Nfp,8) = -1;
        end
        
    end
    
    methods(Access = protected)
        
        function fphys = setInitialField( obj )
            fphys = cell(obj.Nmesh, 1);
            fphys{1}(:,:,1) = ones(size(obj.meshUnion(1).x));
            fphys{1}(:,5,1) = 0;
            fphys{1}(:,:,4) = zeros(size(obj.meshUnion(1).x));
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