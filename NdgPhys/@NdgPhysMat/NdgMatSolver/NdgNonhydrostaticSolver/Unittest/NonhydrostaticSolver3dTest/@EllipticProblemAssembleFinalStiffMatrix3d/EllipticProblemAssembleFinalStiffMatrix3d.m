classdef EllipticProblemAssembleFinalStiffMatrix3d < SWEBarotropic3d
    %ELLIPTICPROBLEMASSEMBLEFINALSTIFFMATRIX3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        
        StiffMatrix
        
        ChLength = 20
        
        ChWidth = 20
        
        Depth = 3.2
        
        SurfaceBoundaryEdgeType = 'Dirichlet'
        
    end
    
    properties(Constant)
        % This parameter is negative such that all randomly produced water
        % depth is larger than the critical depth
        hcrit = -0.01
    end
    
    methods
        function obj = EllipticProblemAssembleFinalStiffMatrix3d(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.varFieldIndex = [1 2 11];
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
        end
        
        function EllipticProblemSolve(obj)
            PHPX = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            PHPY = rand(obj.meshUnion.cell.Np, obj.meshUnion.K);
            
            obj.StiffMatrix = rand(1) * obj.NonhydrostaticSolver.SPNPX + rand(1) * obj.NonhydrostaticSolver.SPNPY;
            FinalStiffMatrix = obj.matAssembleFinalGlobalStiffMatrix( PHPX, PHPY);
            
            obj.NonhydrostaticSolver.TestAssembleFinalGlobalStiffMatrix( obj, obj.fphys, obj.StiffMatrix, PHPX, PHPY);
            disp('======================For StiffMatrix=================================')
            disp('The maximum value is')
            disp(max(max(FinalStiffMatrix)))
            disp('The maximum difference is:')
            disp(max(max(FinalStiffMatrix - obj.NonhydrostaticSolver.GlobalStiffMatrix)))
            disp('The minimum difference is:')
            disp(min(min(FinalStiffMatrix - obj.NonhydrostaticSolver.GlobalStiffMatrix)))
            disp('=====================End StiffMatrix==================================')
        end
    end
    
    methods(Access = protected)
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = rand( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                fphys2d{m}(:,:,1) = obj.Depth*ones( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K) + rand(obj.mesh2d(m).cell.Np, obj.mesh2d(m).K);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 1.75;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('AdvDiffVerticalDiffusionType') = enumVerticalDiffusion.None;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('AdvDiffHorizontalDiffusionType') = enumHorizontalDiffusion.Constant;
            option('AdvDiffConstantHorizontalDiffusionValue') = 1;
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [ 0, obj.ChWidth ], M, M, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

