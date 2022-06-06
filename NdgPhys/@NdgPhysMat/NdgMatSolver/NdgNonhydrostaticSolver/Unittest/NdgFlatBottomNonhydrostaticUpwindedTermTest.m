classdef NdgFlatBottomNonhydrostaticUpwindedTermTest < NdgFlatBottomWetDryTestWithWallBoundary
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgFlatBottomNonhydrostaticUpwindedTermTest(N, cellType)
            obj = obj@NdgFlatBottomWetDryTestWithWallBoundary(N, cellType);
        end
        
        function  [Exactfhx, Exacthux] = GetExactInnerEdgeFluxTermX(obj)
            Exacthux = 2.5*ones(size(obj.meshUnion(1).InnerEdge.nx));
            Exacthux(:,[3 4 6 8 9 11]) = 0;
            Exactfhx = zeros(size(Exacthux));
            %upwinded flux term
            Exactfhx(:,1) = 1;Exactfhx(:,2) = 2; 
            Exactfhx(:,5) = 4;Exactfhx(:,7) = 5;
            Exactfhx(:,10) = 7;Exactfhx(:,12) = 8;  
        end
        
        function [Exactfhx, Exacthux] = GetExactBoundaryEdgeFluxTermX(obj)
            Exacthux = 2.5*ones(size(obj.meshUnion(1).BoundaryEdge.nx));
            Exacthux(:,[1 3 4]) = 0;Exacthux(:,[10 11 12]) = 0;
            Exacthux(:,[2 6 8]) = -1 * Exacthux(:,[2 6 8]);
            Exactfhx = zeros(size(Exacthux));
            Exactfhx(:,2) = -1; Exactfhx(:,6) = -4;Exactfhx(:,8) = -7;
            Exactfhx(:,5) = 3; Exactfhx(:,7) = 6;Exactfhx(:,9) = 9;
        end
        
        
    end
    
    methods(Access = protected)
        
        function fphys = setInitialField( obj )
            fphys = cell(obj.Nmesh, 1);
            fphys{1}(:,:,4) = zeros(size(obj.meshUnion(1).x));
            for i = 1:obj.meshUnion(1).K
                fphys{1}(:,i,1) = i;
            end
            fphys{1}(:,:,2) = 2.5;
            fphys{1}(:,:,3) = 2.5;
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

