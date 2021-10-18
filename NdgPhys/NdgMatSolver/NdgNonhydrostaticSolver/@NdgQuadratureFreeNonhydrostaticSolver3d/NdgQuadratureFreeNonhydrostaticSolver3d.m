classdef NdgQuadratureFreeNonhydrostaticSolver3d < handle
    %NDGQUADRATUREFREENONHYDROSTATICSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        % Second order partial derivative about nonhydrostatic pressure in direction x, $\frac{\partial^2 p}{\partial x^2}$
        SPNPX
        % Second order partial derivative about nonhydrostatic pressure in direction y, $\frac{\partial^2 p}{\partial y^2}$
        SPNPY
        % First order partial derivative about nonhydrostatic pressure in vertical direction, $\frac{\partial p}{\partial \sigma}$
        PNPS
        % First order partial derivative about nonhydrostatic pressure in horizontal direction, $\frac{\partial p}{\partial x}$
        PNPX
        % First order partial derivative about nonhydrostatic pressure in horizontal direction, $\frac{\partial p}{\partial y}$
        PNPY
        % Second order partial derivative about nonhydrostatic pressure in vertical direction, $\frac{\partial^2 p}{\partial \sigma^2}$
        SPNPS
        % Mixed Second order partial derivative about nonhydrostatic pressure in direction x, $\frac{\partial}{\partial x}\left(\frac{\partial p}{\partial \sigma}\right)$
        MSPNPX
        % Mixed Second order partial derivative about nonhydrostatic pressure in direction x, $\frac{\partial}{\partial y}\left(\frac{\partial p}{\partial \sigma}\right)$
        MSPNPY
    end
    
    properties
        % First order partial derivative about $\sigma$ direction x, $\frac{\partial \sigma}{\partial x}$
        PSPX
        % First order partial derivative about $\sigma$ direction y, $\frac{\partial \sigma}{\partial y}$
        PSPY
        % Second order partial derivative about $\sigma$ direction x, $\frac{\partial^2 \sigma}{\partial x^2}$
        SPSPX
        % Second order partial derivative about $\sigma$ direction y, $\frac{\partial^2 \sigma}{\partial y^2}$
        SPSPY
        % Square of the first order partial derivative about $\sigma$ direction x, $\left (\frac{\partial \sigma}{\partial x}\right)^2$
        SQPSPX
        % Square of the first order partial derivative about $\sigma$ direction y, $\left (\frac{\partial \sigma}{\partial y}\right)^2$
        SQPSPY
        % Partial derivative of u with respect to x, $\frac{\partial u}{\partial x}$
        PUPX
        % Partial derivative of u with respect to y, $\frac{\partial u}{\partial y}$
        PUPY
        % Partial derivative of v with respect to x, $\frac{\partial v}{\partial x}$
        PVPX
        % Partial derivative of v with respect to y, $\frac{\partial v}{\partial y}$
        PVPY
        % Partial derivative of u with respect to $\sigma$, $\frac{\partial u}{\partial \sigma}$
        PUPS
        % Partial derivative of v with respect to $\sigma$, $\frac{\partial v}{\partial \sigma}$
        PVPS
        % Partial derivative of w with respect to $\sigma$, $\frac{\partial w}{\partial \sigma}$
        PWPS
        % Partial derivative of H with respect to x, $\frac{\partial H}{\partial x}$
        PHPX
        % Partial derivative of H with respect to y, $\frac{\partial H}{\partial y}$
        PHPY
    end
    
    properties
        GlobalStiffMatrix
        
        BoundNonhydroPressure
        
        BoundNonhydroGrad
        
        NonhydroRHS
        
        Wold
        
        Wnew
        
        Uold
        
        Unew
        
        Vold
        
        Vnew
        
    end
    
    properties
        varIndex
        
        rhsIndex
        
        rho
    end
    
    properties
        mesh
        cell
        InnerEdge
        BoundaryEdge
        BottomBoundaryEdge
        SurfaceBoundaryEdge
        BottomEdge
        mesh2d
        InnerEdge2d
        BoundaryEdge2d
        cell2d
    end
    
    properties
        SortedIEnx
        SortedIEny
        SortedIEGlobalFace
        SortedIEAdjEle
        SortedIEReverseFlag
        SortedIEInternalFace
        UniEleNumber
        UniEle
    end
    
    properties
        PARDISO_INITIALIZED = 0
        PARDISO_INFO = []
    end
    
    methods
        function obj = NdgQuadratureFreeNonhydrostaticSolver3d( PhysClass, mesh )
            obj.matSetInitializeCharacteristicMatrix( PhysClass.SurfaceBoundaryEdgeType, mesh );
            warning('off');
            obj.InnerEdge = struct(mesh.InnerEdge);
            obj.BoundaryEdge = struct(mesh.BoundaryEdge);
            obj.BottomBoundaryEdge = struct(mesh.BottomBoundaryEdge);
            obj.SurfaceBoundaryEdge = struct(mesh.SurfaceBoundaryEdge);
            obj.BottomEdge = struct(mesh.BottomEdge);
            obj.mesh = struct(mesh);
            obj.cell = struct(mesh.cell);
            obj.mesh2d = struct(mesh.mesh2d);
            obj.InnerEdge2d = struct(mesh.mesh2d.InnerEdge);
            obj.BoundaryEdge2d = struct(mesh.mesh2d.BoundaryEdge);
            obj.cell2d = struct(mesh.mesh2d.cell);
            warning('on');
            obj.varIndex = zeros(6,1);
            obj.rho = 1000;
            obj.BoundNonhydroPressure = zeros(obj.BoundaryEdge.Nfp, obj.BoundaryEdge.Ne);
            obj.BoundNonhydroGrad = zeros(obj.BoundaryEdge.Nfp, obj.BoundaryEdge.Ne);
            for i = 1:PhysClass.Nfield
                if (strcmp(PhysClass.fieldName3d{i},'hu'))
                    obj.varIndex(1) = i;
                elseif (strcmp(PhysClass.fieldName3d{i},'hv'))
                    obj.varIndex(2) = i;
                elseif (strcmp(PhysClass.fieldName3d{i},'hw'))
                    obj.varIndex(3) = i;
                elseif (strcmp(PhysClass.fieldName3d{i},'h'))
                    obj.varIndex(4) = i;
                elseif (strcmp(PhysClass.fieldName3d{i},'zx'))
                    obj.varIndex(5) = i;
                elseif (strcmp(PhysClass.fieldName3d{i},'zy'))
                    obj.varIndex(6) = i;
                elseif (strcmp(PhysClass.fieldName3d{i},'z'))
                    obj.varIndex(7) = i;
                end
            end
            obj.rhsIndex = zeros(1,3);
            for i = 1:3
                Index =(find(obj.varIndex(i) == PhysClass.varFieldIndex));
                obj.rhsIndex(i) = Index;
            end
            
            [obj.UniEleNumber, obj.UniEle, obj.SortedIEnx, obj.SortedIEny, obj.SortedIEGlobalFace, ...
                obj.SortedIEAdjEle, obj.SortedIEReverseFlag, obj.SortedIEInternalFace ] = mxGetInnerEdgeTopologyRelation(...
                obj.InnerEdge.FToE, obj.InnerEdge.FToF, obj.InnerEdge.nx, obj.InnerEdge.ny, obj.InnerEdge.Ne, obj.cell.Nface, ...
                obj.InnerEdge.Nfp, obj.mesh.EToE, obj.mesh.K);
            
            obj.matClearGlobalMemory( );
        end
        
        function fphys = NdgConservativeNonhydrostaticUpdata(obj, physClass, fphys, fphys2d, deltatime)
            
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew, obj.Unew, obj.Vnew, obj.PHPX, obj.PHPY ] = ...
                mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            
            edge = physClass.meshUnion.BottomBoundaryEdge;
            [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            obj.Wnew = fm(:,:,obj.varIndex(3))./fm(:,:,obj.varIndex(4));
            
            %             obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrix( obj.PNPS, obj.SPNPX, obj.SPNPY, obj.SPNPS, ...
            %                 obj.MSPNPX, obj.MSPNPY, obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY,...
            %                 fphys{1}(:,:,obj.varIndex(4)), physClass.hcrit, obj.BoundaryEdge2d, ...
            %                 obj.mesh, obj.cell, obj.BoundaryEdge, int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrixNew(obj.SPNPX, obj.SPNPY, obj.PSPX, obj.PSPY, obj.SQPSPX, ...
                obj.SQPSPY, physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.InnerEdge, obj.BottomEdge,...
                obj.SurfaceBoundaryEdge, obj.cell2d.M, obj.mesh2d.J, obj.mesh2d.K, physClass.SurfaceBoundaryEdgeType,...
                obj.SortedIEnx, obj.SortedIEny, obj.SortedIEGlobalFace, obj.SortedIEAdjEle, obj.SortedIEReverseFlag,...
                obj.SortedIEInternalFace);
            
            obj.NonhydroRHS = mxAssembleNonhydroRHS(obj.PUPX, obj.PUPS, obj.PVPY, obj.PVPS, obj.PWPS, obj.PSPX, ...
                obj.PSPY, fphys{1}(:,:,obj.varIndex(4)), deltatime, obj.rho, physClass.hcrit, obj.mesh.J, obj.cell.M);
            
            [  obj.GlobalStiffMatrix, obj.NonhydroRHS ] = mxImposeBoundaryCondition( obj.GlobalStiffMatrix, obj.PSPX, obj.PSPY, ...
                physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.BottomBoundaryEdge,...
                obj.BoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype), obj.mesh2d, obj.cell2d, obj.InnerEdge2d,...
                obj.BoundaryEdge2d, obj.Wold, obj.Wnew, deltatime, obj.rho, fphys{1}(:,:,obj.varIndex(1)), ...
                fphys{1}(:,:,obj.varIndex(2)), obj.NonhydroRHS, obj.PWPS,  obj.BoundNonhydroPressure, obj.Unew, obj.Uold, obj.Vnew, obj.Vold,...
                obj.PUPX, obj.PUPY, obj.PUPS, obj.PVPX, obj.PVPY, obj.PVPS, obj.PHPX, obj.PHPY, physClass.gra, fphys{1}(:,:,obj.varIndex(5)), ...
                fphys{1}(:,:,obj.varIndex(6)), obj.BoundNonhydroGrad);
            
            obj.GlobalStiffMatrix = mxAssembleFinalGlobalStiffMatrix(obj.cell.Np, obj.mesh.K, physClass.hcrit, obj.mesh.EToE, obj.cell.Nface,...
                obj.cell.M, obj.mesh.J, obj.GlobalStiffMatrix, obj.PNPX, obj.PNPY, obj.PNPS, fphys{1}(:,:,obj.varIndex(4)), ...
                obj.PHPX, obj.PHPY, fphys{1}(:,:,obj.varIndex(5)), fphys{1}(:,:,obj.varIndex(6)), obj.mesh.z, obj.UniEleNumber, obj.UniEle);
            
            %             obj.GlobalStiffMatrix = mxAssemblePositiveDefiniteStiffMatrix( obj.GlobalStiffMatrix );
            %
            % %             ittol = 1e-8; maxit = 10000; Droptol = 1e-4;
            % %             Cinc = ichol(obj.GlobalStiffMatrix, struct('droptol', Droptol));
            % %             NonhydroPressure = pcg(obj.GlobalStiffMatrix, -1*obj.NonhydroRHS, ittol, maxit, Cinc', Cinc);
            % %             toc;
            %
            %             NonhydroPressure = obj.GlobalStiffMatrix\(-1*obj.NonhydroRHS);
            
            %             tic;
            %             NonhydroPressure = obj.GlobalStiffMatrix\(obj.NonhydroRHS);
            %             toc;
            %==========================================For unsymmetric matrix======================================================
            %             tic;
            verbose = false;
            if ~obj.PARDISO_INITIALIZED
                
                obj.PARDISO_INFO = pardisoinit(11,0);
                % Analyze the matrix and compute a symbolic factorization.
                obj.PARDISO_INFO = pardisoreorder(obj.GlobalStiffMatrix, obj.PARDISO_INFO,verbose);
                
                obj.PARDISO_INITIALIZED = 1;
            end
            % Compute the numeric factorization.
            obj.PARDISO_INFO = pardisofactor(obj.GlobalStiffMatrix, obj.PARDISO_INFO, false);
            % Compute the solutions X using the symbolic factorization.
            [NonhydroPressure, ~] = pardisosolve(obj.GlobalStiffMatrix, obj.NonhydroRHS, obj.PARDISO_INFO, false);
            
            
            %==========================================For symmetric matrix=======================================================
            %             tic;
            %             obj.GlobalStiffMatrix = mxAssemblePositiveDefiniteStiffMatrix( obj.GlobalStiffMatrix );
            %             verbose = false;
            %             info = pardisoinit(-2,0);
            %             p    = randperm(numel(obj.NonhydroRHS));
            %             info = pardisoreorder(tril(obj.GlobalStiffMatrix),info,verbose,p);
            %             info = pardisofactor(tril(obj.GlobalStiffMatrix),info,verbose);
            %             [SymNonhydroPressure, info] = pardisosolve(tril(obj.GlobalStiffMatrix),obj.NonhydroRHS,info,verbose);
            %             pardisofree(info);
            %             clear info
            %             SymVersion = toc
            
            %==========================================For PETSc usage=============================================================
            %             tic;
            %              [rowptr, colind, val] = crs_matrix(obj.GlobalStiffMatrix);
            %              [TempNonhydroPressure,flag,relres,iter,reshis] = petscSolveCRS(rowptr, colind, val, obj.NonhydroRHS,...
            %                  PETSC_KSPGMRES, 1.e-10, int32(10000000), PETSC_PCJACOBI, 'right');
            %             toc;
            
            
            
            fphys{1}(:,:,obj.varIndex(1:3)) = mxUpdateConservativeFinalVelocity( NonhydroPressure, fphys{1}, obj.varIndex, ...
                obj.rho, deltatime, obj.PSPX, obj.PSPY, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge,...
                obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype));
            
            %             obj.Wold = mxCalculateBottomVerticalVelocity( obj.cell, obj.BottomBoundaryEdge, fphys{1}, obj.varIndex, obj.mesh, physClass.hcrit );
            [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            obj.Uold = fm(:,:,obj.varIndex(1))./fm(:,:,obj.varIndex(4));
            obj.Vold = fm(:,:,obj.varIndex(2))./fm(:,:,obj.varIndex(4));
            obj.Wold = fm(:,:,obj.varIndex(3))./fm(:,:,obj.varIndex(4));
            
            %             obj.NonhydroPressure = obj.NonhydroPressure + reshape(DiffNonhydroPressure, obj.cell.Np, obj.mesh.K);
            
        end
        
        function [ time, fphys2d, fphys3d ] = matHotStart( obj, mesh3d, file2d, file3d, fphys2d, fphys3d)
            ncid2d = netcdf.open(file2d, 'WRITE');
            ncid3d = netcdf.open(file3d, 'WRITE');
            % Time is the first variable in netcdf
            Time = netcdf.getVar(ncid2d, 0);
            time = Time(end);
            % Variable is the second variable in netcdf
            field2d = netcdf.getVar(ncid2d, 1);
            % H, HU and HV
            fphys2d(:,:,1) = field2d(:,:,1, numel(Time));
            fphys2d(:,:,2) = field2d(:,:,2, numel(Time));
            fphys2d(:,:,3) = field2d(:,:,3, numel(Time));
            %Hu, Hv, omega, Hw
            field3d = netcdf.getVar(ncid3d, 1);
            fphys3d(:,:,1) = field3d(:,:,1, numel(Time));
            fphys3d(:,:,2) = field3d(:,:,2, numel(Time));
            fphys3d(:,:,3) = field3d(:,:,3, numel(Time));
            %> H in extended three dimensional fields
            fphys3d(:,:,4) = mesh3d.Extend2dField( fphys2d(:,:,1) );
            fphys3d(:,:,7) = fphys3d(:,:,4) + fphys3d(:,:,6);
            fphys3d(:,:,11) = field3d(:,:,4, numel(Time));
            netcdf.close(ncid2d);
            netcdf.close(ncid3d);
        end
        
        function matCalculateBottomVerticalVelocity( obj, physClass, fphys )
            %             obj.Wold = mxCalculateBottomVerticalVelocity( obj.cell, obj.BottomBoundaryEdge, fphys{1}, obj.varIndex, obj.mesh, physClass.hcrit );
            edge = physClass.meshUnion.BottomBoundaryEdge;
            [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            obj.Uold = fm(:,:,obj.varIndex(1))./fm(:,:,obj.varIndex(4));
            obj.Vold = fm(:,:,obj.varIndex(2))./fm(:,:,obj.varIndex(4));
            obj.Wold = fm(:,:,obj.varIndex(3))./fm(:,:,obj.varIndex(4));
            
        end
        
        function matCalculateNonhydroRHS(obj, physClass, fphys, fphys2d )
            
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew, obj.Unew, obj.Vnew, obj.PHPX, obj.PHPY ] = ...
                mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            physClass.frhs{1}(:,:,obj.rhsIndex) = mxCalculateNonhydroRHS( obj.NonhydroPressure, fphys{1}, obj.rhsIndex, ...
                obj.varIndex, obj.rho, obj.PSPX, obj.PSPY, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge,...
                obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, physClass.frhs{1}, int8(physClass.meshUnion.BoundaryEdge.ftype));
            
        end
        
        function fphys = matUpdataVerticalVelocity( obj, physClass, fphys, fphys2d )
            
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS ] = mxCalculatePartialDerivativeUpdated( physClass.hcrit,...
                obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            %             Rhs = -fphys{1}(:,:,obj.varIndex(4)) .* ( obj.PUPX + obj.PVPY + obj.PUPS .* obj.PSPX + obj.PVPS .* obj.PSPY );
            %
            %             fphys{1}(:,:,obj.varIndex(3)) = fphys{1}(:,:,obj.varIndex(4)) .* physClass.meshUnion(1).VerticalIntegralField(Rhs);
            
            fphys{1}(:,:,obj.varIndex(3)) = mxUpdateVerticalFinalVelocity( fphys{1}, obj.varIndex, physClass.hcrit, obj.PSPX,...
                obj.PSPY, obj.PUPX, obj.PVPY, obj.PUPS, obj.PVPS, physClass.VerticalVelocitySolver.RHSCoeMatrix{1},...
                physClass.VerticalVelocitySolver.VertCoeMatrix{1}, obj.mesh, obj.cell, obj.BottomBoundaryEdge);
        end
        
        function matClearGlobalMemory( obj )
            clear mxCalculatePartialDerivativeUpdated;
            clear mxAssembleGlobalStiffMatrixNew;
            clear mxImposeBoundaryCondition;
            clear mxCalculateBottomVerticalVelocity;
            if obj.PARDISO_INITIALIZED
                obj.PARDISO_INITIALIZED = 0;
                pardisofree(obj.PARDISO_INFO);
                clear obj.PARDISO_INFO
            end
        end
        
        function TestPartialDerivativeCalculation(obj, physClass, fphys, fphys2d)
            
            [ obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, ~ ] = mxCalculatePartialDerivativeUpdated( physClass.hcrit,...
                obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
        end
        
        function RHS = TestAssembleRHS( obj, physClass, fphys, fphys2d, deltatime )
            [ obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS ] = mxCalculatePartialDerivative( physClass.hcrit,...
                obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            RHS = mxAssembleNonhydroRHS(obj.PUPX, obj.PUPS, obj.PVPY, obj.PVPS, obj.PWPS, obj.PSPX, ...
                obj.PSPY, fphys{1}(:,:,obj.varIndex(4)), deltatime, obj.rho, physClass.hcrit);
        end
        
        function Matrix = TestAssembleGlobalStiffMatrix( obj, physClass, fphys, fphys2d )
            
            [ obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS ] = mxCalculatePartialDerivative( physClass.hcrit,...
                obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            Matrix = mxAssembleGlobalStiffMatrix( obj.PNPS, obj.SPNPX, obj.SPNPY, obj.SPNPS, ...
                obj.MSPNPX, obj.MSPNPY, obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY,...
                fphys{1}(:,:,obj.varIndex(4)), physClass.hcrit, obj.BoundaryEdge2d, ...
                obj.mesh, obj.cell, obj.BoundaryEdge, int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
        end
        
        function fphys = TestUpdateConservativeFinalVelocity( obj, physClass, NonhydroPressure, fphys, fphys2d, deltatime)
            
            [ obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.PUVPXY ] = mxCalculatePartialDerivative( physClass.hcrit,...
                obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            fphys{1}(:,:,obj.varIndex(1:3)) = mxUpdateConservativeFinalVelocity( NonhydroPressure, fphys{1}, obj.varIndex, ...
                obj.rho, deltatime, obj.PSPX, obj.PSPY, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge,...
                obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype));
        end
        
        function SimulatedSolution = TestNewFormGlobalStiffMatrix( obj, physClass, fphys )
            %Checked
            obj.PSPX = physClass.K13;
            obj.PSPY = physClass.K23;
            obj.SQPSPX = obj.PSPX .* obj.PSPX;
            obj.SQPSPY = obj.PSPY .* obj.PSPY;
            
            obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrixNew(obj.SPNPX, obj.SPNPY, obj.PSPX, obj.PSPY, obj.SQPSPX, ...
                obj.SQPSPY, physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.InnerEdge, obj.BottomEdge,...
                obj.SurfaceBoundaryEdge, obj.cell2d.M, ...
                obj.mesh2d.J, obj.mesh2d.K, physClass.SurfaceBoundaryEdgeType);
            SimulatedSolution = reshape(obj.GlobalStiffMatrix\physClass.RHS(1:physClass.meshUnion.cell.Np*physClass.meshUnion.K)', physClass.meshUnion.cell.Np, physClass.meshUnion.K);
        end
        
        function [ StiffMatrix, RHS, Wx, Wy, PNPS, PNPX, PNPY ] = TestImposeBoundaryCondition( obj, physClass, fphys, StiffMatrix, RHS, PSPX, PSPY, Wold, Wnew, deltatime, ...
                PWPS, Unew, Uold, Vnew, Vold, PUPX, PUPY, PUPS, PVPX, PVPY, PVPS, PHPX, PHPY)
            % Checked
            [ StiffMatrix, RHS, Wx, Wy, PNPS, PNPX, PNPY  ] = mxImposeBoundaryCondition( StiffMatrix, PSPX, PSPY, ...
                physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.BottomBoundaryEdge,...
                obj.BoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype), obj.mesh2d, obj.cell2d, obj.InnerEdge2d,...
                obj.BoundaryEdge2d, Wold, Wnew, deltatime, obj.rho, fphys{1}(:,:,obj.varIndex(1)), ...
                fphys{1}(:,:,obj.varIndex(2)), RHS, PWPS, obj.BoundNonhydroPressure, Unew, Uold, Vnew, Vold,...
                PUPX, PUPY, PUPS, PVPX, PVPY, PVPS, PHPX, PHPY, physClass.gra, fphys{1}(:,:,obj.varIndex(5)), ...
                fphys{1}(:,:,obj.varIndex(6)));
        end
        
        function TestAssembleFinalGlobalStiffMatrix( obj, physClass, fphys, StiffMatrix, PHPX, PHPY )
            % Checked
            obj.GlobalStiffMatrix = mxAssembleFinalGlobalStiffMatrix(obj.cell.Np, obj.mesh.K, physClass.hcrit, obj.mesh.EToE, obj.cell.Nface,...
                obj.cell.M, obj.mesh.J, StiffMatrix, obj.PNPX, obj.PNPY, obj.PNPS, fphys{1}(:,:,obj.varIndex(4)), ...
                PHPX, PHPY, fphys{1}(:,:,obj.varIndex(5)), fphys{1}(:,:,obj.varIndex(6)), obj.mesh.z);
        end
        
    end
    
end

