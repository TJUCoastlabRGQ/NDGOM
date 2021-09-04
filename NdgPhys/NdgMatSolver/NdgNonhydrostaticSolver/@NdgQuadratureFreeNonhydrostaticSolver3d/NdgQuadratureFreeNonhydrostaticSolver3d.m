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
        % Partial derivative of v with respect to x, $\frac{\partial v}{\partial y}$
        PVPY
        % Partial derivative of u with respect to $\sigma$, $\frac{\partial u}{\partial \sigma}$
        PUPS
        % Partial derivative of v with respect to $\sigma$, $\frac{\partial v}{\partial \sigma}$
        PVPS
        % Partial derivative of w with respect to $\sigma$, $\frac{\partial w}{\partial \sigma}$
        PWPS
    end
    
    properties
        GlobalStiffMatrix
        
        BoundNonhydroPressure
        
        NonhydroRHS
        
        Wold
        
        Wnew
    end
    
    properties
        varIndex
        
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
            obj.matClearGlobalMemory( );
        end
        
        function fphys = NdgConservativeNonhydrostaticUpdata(obj, physClass, fphys, fphys2d, deltatime)
            
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew ] = mxCalculatePartialDerivativeUpdated( physClass.hcrit,...
                obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            %             edge = physClass.meshUnion.BottomBoundaryEdge;
            %             [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            %             obj.Wnew = fm(:,:,obj.varIndex(3))./fm(:,:,obj.varIndex(4));
            
            %             obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrix( obj.PNPS, obj.SPNPX, obj.SPNPY, obj.SPNPS, ...
            %                 obj.MSPNPX, obj.MSPNPY, obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY,...
            %                 fphys{1}(:,:,obj.varIndex(4)), physClass.hcrit, obj.BoundaryEdge2d, ...
            %                 obj.mesh, obj.cell, obj.BoundaryEdge, int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
            
            obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrixNew(obj.SPNPX, obj.SPNPY, obj.PSPX, obj.PSPY, obj.SQPSPX, ...
                obj.SQPSPY, physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.InnerEdge, obj.BottomEdge,...
                obj.SurfaceBoundaryEdge, obj.cell2d.M, ...
                obj.mesh2d.J, obj.mesh2d.K, physClass.SurfaceBoundaryEdgeType);
            
            obj.NonhydroRHS = mxAssembleNonhydroRHS(obj.PUPX, obj.PUPS, obj.PVPY, obj.PVPS, obj.PWPS, obj.PSPX, ...
                obj.PSPY, fphys{1}(:,:,obj.varIndex(4)), deltatime, obj.rho, physClass.hcrit, obj.mesh.J, obj.cell.M);
            
            [  obj.GlobalStiffMatrix, obj.NonhydroRHS ] = mxImposeBoundaryCondition( obj.GlobalStiffMatrix, obj.PSPX, obj.PSPY, ...
                physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.BottomBoundaryEdge,...
                obj.BoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype), obj.mesh2d, obj.cell2d, obj.InnerEdge2d,...
                obj.BoundaryEdge2d, obj.Wold, obj.Wnew, deltatime, obj.rho, fphys{1}(:,:,obj.varIndex(1)), ...
                fphys{1}(:,:,obj.varIndex(2)), obj.NonhydroRHS, obj.PWPS,  obj.BoundNonhydroPressure);
            
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
            tic;
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
            toc;

%             pardisofree(PARDISO_INFO);
%             clear PARDISO_INFO
            
            
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
            
            obj.Wold = mxCalculateBottomVerticalVelocity( obj.cell, obj.BottomBoundaryEdge, fphys{1}, obj.varIndex, obj.mesh, physClass.hcrit );
            
            %             fphys{1}(1:4,:,obj.varIndex(3)) = 0;
        end
        
        function matCalculateBottomVerticalVelocity( obj, physClass, fphys )
            obj.Wold = mxCalculateBottomVerticalVelocity( obj.cell, obj.BottomBoundaryEdge, fphys{1}, obj.varIndex, obj.mesh, physClass.hcrit );
            
            %             edge = physClass.meshUnion.BottomBoundaryEdge;
            %             [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            %             obj.Wold = fm(:,:,obj.varIndex(3))./fm(:,:,obj.varIndex(4));
            
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
        
    end
    
end

