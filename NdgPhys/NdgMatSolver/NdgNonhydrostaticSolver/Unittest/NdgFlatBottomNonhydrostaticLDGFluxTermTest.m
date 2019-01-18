classdef NdgFlatBottomNonhydrostaticLDGFluxTermTest < NdgFlatBottomWetDryTestWithWallBoundary
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgFlatBottomNonhydrostaticLDGFluxTermTest(N, cellType)
            obj = obj@NdgFlatBottomWetDryTestWithWallBoundary(N, cellType);
        end
        
        function [ExactScalerJumpX, ExactScalerJumpY] = getExactHuFieldInnerEdgeScalerJump(obj)
             ExactScalerJumpX = zeros(size(obj.meshUnion(1).InnerEdge.nx)); 
             ExactScalerJumpY = zeros(size(obj.meshUnion(1).InnerEdge.nx)); 
             ExactScalerJumpX(:,1) = -3;ExactScalerJumpX(:,2) = -5;
             ExactScalerJumpX(:,5) = -9;ExactScalerJumpX(:,7) = -11;
             ExactScalerJumpX(:,10) = -15;ExactScalerJumpX(:,12) = -17;
             ExactScalerJumpY(:,3) = -15;ExactScalerJumpY(:,4) = -21;
             ExactScalerJumpY(:,6) = -27;ExactScalerJumpY(:,8) = -33;
             ExactScalerJumpY(:,9) = -39;ExactScalerJumpY(:,11) = -45;
        end
        
        function [ExactScalerJumpX, ExactScalerJumpY] = getExactHvFieldInnerEdgeScalerJump(obj)
             ExactScalerJumpX = zeros(size(obj.meshUnion(1).InnerEdge.ny)); 
             ExactScalerJumpY = zeros(size(obj.meshUnion(1).InnerEdge.ny));   
             ExactScalerJumpY(:,3) = -63;ExactScalerJumpY(:,4) = -117;
             ExactScalerJumpY(:,6) = -189;ExactScalerJumpY(:,8) = -279;
             ExactScalerJumpY(:,9) = -387;ExactScalerJumpY(:,11) = -513;
             ExactScalerJumpX(:,1) = -7;ExactScalerJumpX(:,2) = -19;
             ExactScalerJumpX(:,5) = -61;ExactScalerJumpX(:,7) = -91;
             ExactScalerJumpX(:,10) = -169;ExactScalerJumpX(:,12) = -217;             
        end
        
        function ExactVectorJump = getExactInnerEdgeVectorJump(obj)
            ExactVectorJump = zeros(size(obj.meshUnion(1).InnerEdge.ny)); 
            ExactVectorJump(:,1) = -3;ExactVectorJump(:,2) = -5;
            ExactVectorJump(:,5) = -9;ExactVectorJump(:,7) = -11;
            ExactVectorJump(:,10) = -15;ExactVectorJump(:,12) = -17;
            ExactVectorJump(:,3) = -63;ExactVectorJump(:,4) = -117;
            ExactVectorJump(:,6) = -189;ExactVectorJump(:,8) = -279;
            ExactVectorJump(:,9) = -387;ExactVectorJump(:,11) = -513;  
        end
        
        function [ ExactScalerJumpX, ExactScalerJumpY ] = getExactHuFieldBoundaryEdgeScalerJump(obj, ftype)
             ExactScalerJumpX = zeros(size(obj.meshUnion(1).InnerEdge.ny)); 
             ExactScalerJumpY = zeros(size(obj.meshUnion(1).InnerEdge.ny));               
            if ftype == enumNonhydroBoundaryCondition.Zero
                % doing nothing
            elseif ftype == enumNonhydroBoundaryCondition.ZeroGrad
                ExactScalerJumpX(:,2) = -2; ExactScalerJumpY(:,1) = -2;
                ExactScalerJumpX(:,5) = 18; ExactScalerJumpY(:,3) = -8;
                ExactScalerJumpX(:,6) = -32; ExactScalerJumpY(:,4) = -18;
                ExactScalerJumpX(:,7) = 72; ExactScalerJumpY(:,10) = 98;
                ExactScalerJumpX(:,8) = -98; ExactScalerJumpY(:,11) = 128;
                ExactScalerJumpX(:,9) = 162; ExactScalerJumpY(:,12) = 162;
            end
        end
        
        function [ ExactScalerJumpX, ExactScalerJumpY ] = getExactHvFieldBoundaryEdgeScalerJump(obj, ftype)
             ExactScalerJumpX = zeros(size(obj.meshUnion(1).InnerEdge.ny)); 
             ExactScalerJumpY = zeros(size(obj.meshUnion(1).InnerEdge.ny));               
            if ftype == enumNonhydroBoundaryCondition.Zero

            elseif ftype == enumNonhydroBoundaryCondition.ZeroGrad
                ExactScalerJumpX(:,2) = -2;  ExactScalerJumpY(:,1) = -2; 
                ExactScalerJumpX(:,5) = 54;  ExactScalerJumpY(:,3) = -16; 
                ExactScalerJumpX(:,6) = -128; ExactScalerJumpY(:,4) = -54;
                ExactScalerJumpX(:,7) = 432;  ExactScalerJumpY(:,10) = 686;
                ExactScalerJumpX(:,8) = -686; ExactScalerJumpY(:,11) = 1024;
                ExactScalerJumpX(:,9) = 1458; ExactScalerJumpY(:,12) = 1458;
            end            
        end
        
        function ExactVectorJump = getExactBoundaryEdgeVectorJump(obj, ftype)
            ExactVectorJump = zeros(size(obj.meshUnion(1).InnerEdge.ny));  
            if ftype == enumNonhydroBoundaryCondition.Zero
               
            elseif ftype == enumNonhydroBoundaryCondition.ZeroGrad
                ExactVectorJump(:,1) = -2;ExactVectorJump(:,2) = -2;
                ExactVectorJump(:,3) = -16;ExactVectorJump(:,4) = -54;
                ExactVectorJump(:,5) = 18;ExactVectorJump(:,6) = -32;
                ExactVectorJump(:,7) = 72;ExactVectorJump(:,8) = -98;
                ExactVectorJump(:,9) = 162;ExactVectorJump(:,10) = 686;
                ExactVectorJump(:,11) = 1024;ExactVectorJump(:,12) = 1458;
            end                 
        end
        
    end
    
    methods(Access = protected)
        
        function fphys = setInitialField( obj )
            fphys = cell(obj.Nmesh, 1);
            fphys{1}(:,:,4) = zeros(size(obj.meshUnion(1).x));
            for i = 1:obj.meshUnion(1).K
                fphys{1}(:,i,1) = i;
                fphys{1}(:,i,2) = i.^2;
                fphys{1}(:,i,3) = i.^3;
            end
            
            
        end
        
        function [ mesh ] = makeUniformMesh(obj, N, type)
            bctype = [...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, [-3, 3], [-3, 3], 3, 3, bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, [-3, 3], [-3, 3], 3, 3, bctype);
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func
    end
    
end

