classdef ModeSplitStandingWaveInAClosedChannel2d < SWEPreBlanaced2d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 1e-3
    end
    
    
    methods
        function obj = ModeSplitStandingWaveInAClosedChannel2d( mesh )
                 obj = obj@SWEPreBlanaced2d();
                 obj.initPhysFromOptions( mesh );
        end
        
%         AnalysisResult2d( obj );
%         AnalysisResult3d( obj );
        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys2d{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys2d{m}(:,:,1) = ones( mesh.cell.Np, mesh.K );
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 500;                  
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
        end
        
    end
end

