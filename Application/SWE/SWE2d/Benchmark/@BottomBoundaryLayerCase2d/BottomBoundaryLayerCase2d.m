classdef BottomBoundaryLayerCase2d < SWEConventional2d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 10000;
        %> channel depth
        H0 = 15;
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 86400;
        % to be corrected
    end
    
    properties
        dt
    end
    
    methods
        function obj = BottomBoundaryLayerCase2d(N, deltax, cellType)
            obj = obj@SWEConventional2d();
            obj.hmin = 1e-3;
            [ mesh ] = obj.makeUniformMesh(N, deltax, cellType);
            obj.initPhysFromOptions( mesh );
            obj.outputFieldOrder = [1, 2, 3, 6];   
        end
        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [ fphys] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                % bottom elevation
                fphys{m}(:, :, 4) = -obj.H0;                
                %water depth
                fphys{m}(:,:,1) = -obj.meshUnion(1).x *10^(-5) - fphys{m}(:, :, 4);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 86400;
            outputIntervalNum = 2500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK33;
            option('outputNcfileNum') = 500;            
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
        end
        
        function [ mesh ] = makeUniformMesh(obj, N, deltax, type)
            bctype = [...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, [-obj.ChLength/2, obj.ChLength/2], 0.1*[-obj.ChLength/2, obj.ChLength/2], obj.ChLength/deltax, 0.1*obj.ChLength/deltax, bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, [-obj.ChLength/2, obj.ChLength/2], 0.1*[-obj.ChLength/2, obj.ChLength/2], obj.ChLength/deltax, 0.1*obj.ChLength/deltax, bctype);% 20/0.1 22/0.05  %4/0.025, 1/0.0125,
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func  
        
    end
end


