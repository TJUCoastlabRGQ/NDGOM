classdef NdgFlatBottomForSecondOrderTerm  < NdgNonhydrostaticAbstractTest
%<@brief function used to test the calculation for the second order term contained in the left hand side of the nonhydrostatic model
    
    properties
    end
    
    methods
        function obj = NdgFlatBottomForSecondOrderTerm(N, cellType)
            obj = obj@NdgNonhydrostaticAbstractTest(N, cellType);
        end      
        
        function [qx, qy, q2x, q2y] = getExactSecondOrderTerm(obj)
            mesh = obj.meshUnion(1);
            qx = 2 * ( mesh.x + 5 );qy = zeros(size(mesh.x));
            q2x = 2 * ones(size(mesh.x)); q2y = zeros(size(mesh.x)); 
        end
        
    end

    methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = -10;
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
                fphys{m}(:,:,1) =  ( mesh.x + 5 ) .^ 2;
%                 fphys{m}(:,:,1) =  obj.A * cos(2*pi*mesh.x/obj.Lambda) * cos( 2*pi * sqrt(obj.gra*obj.d)/obj.Lambda*0) - fphys{m}(:,:,4);                
            end
        end
        
        function [ mesh ] = makeUniformMesh(obj, N, type)
            bctype = [...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, [-2, 2], [-1, 1], 2, 1, bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, [0, 40], [0, 2], 20, 2, bctype);
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func        
    end    
end

