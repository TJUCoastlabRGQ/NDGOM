classdef SecondOrderOperatorTest < SWEAbstract3d
    %SECONDORDEROPERATORTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ChLength = 200
        ChWidth = 200
        miu = 0.01
        Nfield2d = 1
        Nvar2d = 1
        varFieldIndex2d = 1
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    properties
        outputFieldOrder2d = []
        outputFieldOrder =  1
        Nfield = 1
        Nvar = 1
        varFieldIndex = 1
        fieldName3d = {'C'};
    end
    
    properties
        %         mesh2d
        %         mesh3d
    end
    
    methods
        function obj = SecondOrderOperatorTest(N, Nz, M, Mz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.meshUnion = obj.mesh3d;
            obj.Nmesh = 1;
            %             obj.outputFieldOrder2d = [];
            obj.fphys = obj.setInitialField;
            obj.option = obj.setOption( obj.option );
            obj.outputFile = obj.matInitOutput;
        end
        
        matTimeStepping(obj);
    end
    
    methods( Hidden )
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            %an empty function
        end
    end
    
    methods(Access = protected)
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            fphys{1}(:,:,1) = 1/sqrt(4*0+1)*exp(-(obj.mesh3d(1).z+0.5).^2/obj.miu/(4*0+1));
%             fphys{1}(:,:,1) = 1/obj.miu*exp(-(obj.mesh3d(1).z+0.5).^2);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 2;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('EddyViscosityType') = enumEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantEddyViscosityValue') = 0.01;
        end
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformTriMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

