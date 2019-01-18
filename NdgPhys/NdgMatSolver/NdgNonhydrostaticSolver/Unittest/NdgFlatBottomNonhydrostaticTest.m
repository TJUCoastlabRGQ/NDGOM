classdef NdgFlatBottomNonhydrostaticTest < NdgNonhydrostaticAbstractTest
    %NDGFLATBOTTOMNONHYDROSTATICTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgFlatBottomNonhydrostaticTest(N, cellType)
            obj = obj@NdgNonhydrostaticAbstractTest(N, cellType);
        end  
        
        function ExactEidBoundaryType = getExactEidBoundaryType(obj)
            ExactEidBoundaryType = 1 * ones(size(obj.meshUnion(1).BoundaryEdge.FToN1));
        end
    end
    
   methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell(obj.Nmesh, 1);
            fphys{1}(:,:,4) = zeros(size(obj.meshUnion(1).x));
            fphys{1}(:,:,1) = 4 - fphys{1}(:,:,4);
        end
        
        function [ mesh ] = makeUniformMesh(obj, N, type)
            bctype = [...
                enumBoundaryCondition.ClampedDepth, ...
                enumBoundaryCondition.ClampedDepth, ...
                enumBoundaryCondition.ClampedDepth, ...
                enumBoundaryCondition.ClampedDepth];
            
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

