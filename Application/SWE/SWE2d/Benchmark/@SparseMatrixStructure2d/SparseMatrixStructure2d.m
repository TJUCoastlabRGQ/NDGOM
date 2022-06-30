classdef SparseMatrixStructure2d < SWEPreBlanaced2d
    %SPARSEMATRIXSTRUCTURE2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
   properties
    SMSFile = 'Application\SWE\SWE2d\Benchmark\@SparseMatrixStructure2d\fort.14'
    end
    
    methods
        function obj = SparseMatrixStructure2d( N )
            obj = obj@SWEPreBlanaced2d();
            [ mesh ] = makeSMSMesh(N, obj.SMSFile);
            obj.hmin = 1e-3;      
            obj.initPhysFromOptions( mesh );
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
                fphys{m}(:,:,1) = -bot;
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
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
            option('outputNcfileNum') = 1;
            option('outputType') = enumOutputFile.VTK;
%             option('nonhydrostaticType') = enumNonhydrostaticType.Hydrostatic;
        end
    end
    
end

function [ mesh ] = makeSMSMesh(N, filename)

mesh = makeSMSFileUMeshUnion2d( N, filename );
end% func
