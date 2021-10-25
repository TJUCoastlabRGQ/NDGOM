classdef NonhydroImposeBCsTest3d < SWEBarotropic3d
    %NONHYDROIMPOSEBCSTEST3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        StiffMatrix
        
        RHS
        
        ChLength = 20
        
        ChWidth = 20
        
        Depth = 3.2
        
        SurfaceBoundaryEdgeType = 'Dirichlet'
        
        MatWx
        
        MatWy
        
        MatPNPS
        
        MatPNPX
        
        MatPNPY
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    methods
        function obj = NonhydroImposeBCsTest3d(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.varFieldIndex = [1 2 11];
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
        end
        
        function EllipticProblemSolve(obj)
            
            TempStiffMatrix = abs(rand(obj.meshUnion.cell.Np*obj.meshUnion.K));
            
            TempRHS = rand(obj.meshUnion.cell.Np*obj.meshUnion.K, 1);
            PSPX = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PSPY = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            Wold = rand(obj.meshUnion.mesh2d.cell.Np,obj.meshUnion.mesh2d.K);
            deltatime = 0.05;
            PWPS = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            Uold = rand(obj.meshUnion.mesh2d.cell.Np,obj.meshUnion.mesh2d.K);
%             Unew =
            edge = obj.meshUnion.BottomBoundaryEdge;
            [fm, ~] = edge.matEvaluateSurfValue(obj.fphys);
            Unew = fm(:,:,1)./fm(:,:,4);
            Vold = rand(obj.meshUnion.mesh2d.cell.Np,obj.meshUnion.mesh2d.K);
            Vnew = fm(:,:,2)./fm(:,:,4);
            Wnew = fm(:,:,11)./fm(:,:,4);
            PUPX = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PUPY = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PUPS = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PVPX = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PVPY = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PVPS = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PHPX = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            PHPY = rand(obj.meshUnion.cell.Np,obj.meshUnion.K);
            obj.matAssembleStiffMatrixAndRHS( TempRHS, PSPX, PSPY, Wold, Wnew, deltatime, ...
                PWPS, Unew, Uold, Vnew, Vold, PUPX, PUPY, PUPS, PVPX, PVPY, PVPS, PHPX, PHPY);
            
            [ ~, NewRHS] = obj.NonhydrostaticSolver.TestImposeBoundaryCondition( obj, ...
                obj.fphys, sparse(TempStiffMatrix), TempRHS, PSPX, PSPY, Wold, Wnew, deltatime, ...
                PWPS, Unew, Uold, Vnew, Vold, PUPX, PUPY, PUPS, PVPX, PVPY, PVPS, PHPX, PHPY);
            
            disp('==================================For RHS========================================')
            disp('The maximum difference is:')
            disp(max(max(NewRHS - obj.RHS)));
            disp('The minimum difference is:')
            disp(min(min(NewRHS - obj.RHS)))
            disp('=================================End RHS========================================')
            
%             disp('==================================For Wx========================================')
%             disp('The maximum value is:')
%             disp(max(max(Wx)))
%             disp('The maximum difference is:')
%             disp(max(max(Wx - obj.MatWx)));
%             disp('The minimum difference is:')
%             disp(min(min(Wx - obj.MatWx)))
%             disp('=================================End Wx========================================')
%             
%             disp('==================================For Wy========================================')
%             disp('The maximum value is:')
%             disp(max(max(Wy)))
%             disp('The maximum difference is:')
%             disp(max(max(Wy - obj.MatWy)));
%             disp('The minimum difference is:')
%             disp(min(min(Wy - obj.MatWy)))
%             disp('=================================End Wy========================================')
%             
%             disp('==================================For PNPS========================================')
%             disp('The maximum value is:')
%             disp(max(max(PNPS)))
%             disp('The maximum difference is:')
%             disp(max(max(PNPS - obj.MatPNPS)));
%             disp('The minimum difference is:')
%             disp(min(min(PNPS - obj.MatPNPS)))
%             disp('=================================End PNPS========================================')
%             
%             disp('==================================For PNPX========================================')
%             disp('The maximum value is:')
%             disp(max(max(PNPX)))
%             disp('The maximum difference is:')
%             disp(max(max(PNPX - obj.MatPNPX)));
%             disp('The minimum difference is:')
%             disp(min(min(PNPX - obj.MatPNPX)))
%             disp('=================================End PNPX========================================')
%             
%             disp('==================================For PNPY========================================')
%             disp('The maximum value is:')
%             disp(max(max(PNPY)))
%             disp('The maximum difference is:')
%             disp(max(max(PNPY - obj.MatPNPY)));
%             disp('The minimum difference is:')
%             disp(min(min(PNPY - obj.MatPNPY)))
%             disp('=================================End PNPY========================================')
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