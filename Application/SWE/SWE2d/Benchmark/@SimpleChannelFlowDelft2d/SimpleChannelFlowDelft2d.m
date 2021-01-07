classdef SimpleChannelFlowDelft2d < SWEPreBlanaced2d
    %SIMPLECHANNELFLOWDELFT2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Length = 10000
        
        Width = 2000
    end
    
    methods
        function obj = SimpleChannelFlowDelft2d( N, deltax, cellType )
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;
            [ mesh ] = obj.makeUniformMesh(N, deltax, cellType);
            obj.initPhysFromOptions( mesh );
        end
    end
    
    methods(Access = protected)
        
        %> set initial function
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys{m} = zeros( obj.meshUnion.cell.Np, obj.meshUnion.K, obj.Nfield );
                % bottom elevation
                fphys{m}(:, :, 4) =  -0.0001 * obj.meshUnion.x;
                %water depth
                fphys{m}(:,:,1) = 3.89677 - 0.0001 * (10000-obj.meshUnion.x);
                %                  fphys2d{m}(:,:,1) = 2.89677;
            end
        end
        
        function matUpdateExternalField( obj, time, ~ )
            hu2d = zeros(size(obj.fext{1}(:,:,1)));
            h2d  = zeros(size(obj.fext{1}(:,:,1)));
            Index = ( obj.meshUnion.BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel );
            hu2d(:,Index) = 5;
            obj.fext{1}(:,:,2) = hu2d;
            Index = ( obj.meshUnion.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
            h2d(:,Index) = 3.89677;
            obj.fext{1}(:,:,1) = h2d;
        end
        
        
        function [ mesh ] = makeUniformMesh(obj, N, deltax, type)
            bctype = [ ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.ClampedVel, ...
                enumBoundaryCondition.ClampedDepth ];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, [0, obj.Length], [0, obj.Width], ceil( obj.Length/deltax ), ceil( obj.Width/deltax ), bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, [0, obj.Length], [0, obj.Width], ceil( obj.Length/deltax ), ceil( obj.Width/deltax ), bctype);
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func
        
        function [ option ] = setOption( obj, option )
            ftime = 86400;
            outputIntervalNum = 3000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK33;
            option('outputNcfileNum') = 500;
            option('limiterType') = enumLimiter.Vert;
            option('outputType') = enumOutputFile.VTK;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
        end
    end
    
end

