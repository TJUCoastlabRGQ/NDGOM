classdef TideInSemiClosedChannel < SWEPreBlanaced2d
    %TIDEINSEMICLOSEDCHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        amplitude = 0.25
        Depth = 30
        T = 12*3600
%         omega = 2*pi/12/3600
        ChLength = 24000
        ChWidth = 2000
    end
    
    properties
        Length
        k
    end
    
    methods
        function obj = TideInSemiClosedChannel(N, deltax, cellType)
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;
            [ mesh ] = obj.makeUniformMesh(N, deltax, cellType);
            obj.initPhysFromOptions( mesh );
            obj.WaveCharacterEstimate;
            obj.outputFieldOrder = [1, 2, 3];
        end
        
    end
    
    methods(Access = protected)

        function WaveCharacterEstimate(obj)
            f = @(L) L - obj.gra*(obj.T)^2/(2*pi)*tanh(2*pi/L*obj.Depth);
            obj.Length = fzero(f,[0.01 20000000]);
        end
        
        function matUpdateExternalField( obj, time, ~ )
            
            obj.fext{1}( :, :, 1 ) = obj.amplitude * sin(2*pi/obj.T*time )+obj.Depth;
            
        end
        
        function fphys = setInitialField( obj )
            fphys = cell( 1, 1 );
            mesh = obj.meshUnion(1);
            fphys{1} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
            fphys{1}(:,:,1) = obj.Depth;
            
            fphys{1}(:,:,4) = -fphys{1}(:,:,1);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 21600;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.EXRK33;
            option('outputNcfileNum') = 1;
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
                enumBoundaryCondition.ClampedDepth];
            
            if (type == enumStdCell.Tri)
                mesh = makeUniformTriMesh(N, [0, obj.ChLength], [-obj.ChWidth/2, obj.ChWidth/2], obj.ChLength/deltax, obj.ChWidth/deltax, bctype);
            elseif(type == enumStdCell.Quad)
                mesh = makeUniformQuadMesh(N, [0, obj.ChLength], [-obj.ChWidth/2, obj.ChWidth/2], obj.ChLength/deltax, obj.ChWidth/deltax, bctype);
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func
        
    end
    
end

