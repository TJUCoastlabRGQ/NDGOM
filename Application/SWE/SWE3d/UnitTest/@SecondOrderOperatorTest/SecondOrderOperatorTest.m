classdef SecondOrderOperatorTest < SWEAbstract3d
    %SECONDORDEROPERATORTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ChLength = 20
        ChWidth = 20
        Nfield2d = 1
        Nvar2d = 1
        varFieldIndex2d = 1
    end
    
    properties
        miu = 0.01
        
        Cexact
        
        DiffCexact
        
        DirichExact
        
        NewmannExact
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    properties
        outputFieldOrder2d = []
        outputFieldOrder3d =  1
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
            obj.matGetFunction;
            %             obj.outputFieldOrder2d = [];
            obj.fphys = obj.setInitialField;
            obj.option = obj.setOption( obj.option );
            obj.outputFile3d = obj.matInitOutput(obj.mesh3d, obj.fieldName3d);
            obj.DirichExact = zeros(obj.mesh2d.cell.Np,2);
            obj.NewmannExact = zeros(obj.mesh2d.cell.Np,2);
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
            z = obj.meshUnion.z;
            t = 0;
            fphys{1}(:,:,1) = eval(obj.Cexact).*ones(size(obj.meshUnion.x));
            %             fphys{1}(:,:,1) = 1/obj.miu*exp(-(obj.mesh3d(1).z+0.5).^2);
        end
        
        function matUpdateExternalField( obj, time )
            %Top first column, bottom second column
            t = time;
            z = zeros(obj.mesh2d.cell.Np,1);
            obj.DirichExact(:,1) = eval(obj.Cexact);
            obj.NewmannExact(:,1) = obj.miu * eval(obj.DiffCexact);
            %                         obj.NewmannExact(1) = obj.miu*1;
            z = -1 * ones(obj.mesh2d.cell.Np,1);
            obj.DirichExact(:,2) = eval(obj.Cexact);
            obj.NewmannExact(:,2) = obj.miu * eval(obj.DiffCexact);
        end
        
        function matGetFunction(obj)
            syms z t;
            z0 = -0.5;
            obj.Cexact = 1/sqrt(4*t+1)*exp(-(z-z0)^2/obj.miu/(4*t+1));
            obj.DiffCexact = diff(obj.Cexact, z);
            %                         obj.Cexact = 0*t+0.*x;
            %             obj.Cexact =-1* sin(2*pi*x + 0 * t);
            %             obj.DiffCexact = diff(obj.Cexact, x);
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 400;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantVerticalEddyViscosityValue') = 0.01;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 100;
        end
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

