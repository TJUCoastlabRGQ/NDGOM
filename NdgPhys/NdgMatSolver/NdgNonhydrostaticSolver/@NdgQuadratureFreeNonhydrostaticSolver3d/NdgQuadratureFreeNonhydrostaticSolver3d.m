classdef NdgQuadratureFreeNonhydrostaticSolver3d < handle
    %NDGQUADRATUREFREENONHYDROSTATICSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        % First order partial derivative about nonhydrostatic pressure in vertical direction, $\frac{\partial p}{\partial \sigma}$
        PNPS
        % First order partial derivative about nonhydrostatic pressure in horizontal direction, $\frac{\partial p}{\partial x}$
        PNPX
        % First order partial derivative about nonhydrostatic pressure in horizontal direction, $\frac{\partial p}{\partial y}$
        PNPY
    end
    
    properties
        % First order partial derivative about $\sigma$ direction x, $\frac{\partial \sigma}{\partial x}$
        PSPX
        % First order partial derivative about $\sigma$ direction y, $\frac{\partial \sigma}{\partial y}$
        PSPY
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
        
        BotBEUold
        
        BotBEUnew
        
        BotBEVold
        
        BotBEVnew
        
        BotBEWold
        
        BotBEWnew
        
        BEUold
        
        BEUnew
        
        BEVold
        
        BEVnew
        
        BEWold
        
        BEWnew
        
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
    
    methods
        function obj = NdgQuadratureFreeNonhydrostaticSolver3d( PhysClass, mesh )
            
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
            
            obj.matSetInitializeCharacteristicMatrix( mesh );
            
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
            
            %checked
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew, obj.Unew, obj.Vnew, obj.PHPX, obj.PHPY ] = ...
                mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            
            obj.matGetBottomNewMomentum( physClass, fphys );
            
            % checked
           [ obj.GlobalStiffMatrix, NNZ ] = mxAssembleGlobalStiffMatrixNew(obj.PSPX, obj.PSPY, obj.SQPSPX, ...
                obj.SQPSPY, physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.InnerEdge, obj.BottomEdge,...
                obj.SurfaceBoundaryEdge, obj.cell2d.M, obj.mesh2d.J, obj.mesh2d.K, physClass.SurfaceBoundaryEdgeType,...
                obj.SortedIEnx, obj.SortedIEny, obj.SortedIEGlobalFace, obj.SortedIEAdjEle, obj.SortedIEReverseFlag,...
                obj.SortedIEInternalFace, obj.BoundaryEdge, int8(obj.BoundaryEdge.ftype));
            %checked
            obj.NonhydroRHS = mxAssembleNonhydroRHS(obj.PUPX, obj.PUPS, obj.PVPY, obj.PVPS, obj.PWPS, obj.PSPX, ...
                obj.PSPY, fphys{1}(:,:,obj.varIndex(4)), deltatime, obj.rho, physClass.hcrit, obj.mesh.J, obj.cell.M);
            
            
            %             obj.NonhydroRHS = mxImposeBottomBoundaryCondition( obj.PSPX, obj.PSPY, physClass.hcrit, obj.mesh, obj.cell,...
            %                 obj.BottomBoundaryEdge, obj.BotBEUold, obj.BotBEUnew, obj.BotBEVold, obj.BotBEVnew, obj.BotBEWold, obj.BotBEWnew,...
            %                 deltatime, obj.rho, obj.varIndex, obj.NonhydroRHS, physClass.frhs{1},...
            %                 fphys{1});
            
            %             obj.NonhydroRHS = mxImposeSideBoundaryCondition( obj.PSPX, obj.PSPY, physClass.hcrit, obj.mesh, obj.cell,...
            %                 obj.BoundaryEdge, obj.BEUold, obj.BEUnew, obj.BEVold, obj.BEVnew, obj.BEWold, obj.BEWnew,...
            %                 deltatime, obj.rho, obj.varIndex, obj.NonhydroRHS, physClass.frhs{1},...
            %                 fphys{1}, int8(obj.BoundaryEdge.ftype));
            
            
            % checked
            obj.GlobalStiffMatrix = mxAssembleFinalGlobalStiffMatrix(obj.cell.Np, obj.mesh.K, physClass.hcrit, obj.mesh.EToE, obj.cell.Nface,...
                obj.cell.M, obj.mesh.J, obj.GlobalStiffMatrix, obj.PNPX, obj.PNPY, obj.PNPS, fphys{1}(:,:,obj.varIndex(4)), ...
                obj.PHPX, obj.PHPY, fphys{1}(:,:,obj.varIndex(5)), fphys{1}(:,:,obj.varIndex(6)), obj.mesh.z, obj.UniEleNumber, obj.UniEle);
 
            NonhydroPressure = mxNonhydroSystemSolve( obj.GlobalStiffMatrix, NNZ, obj.cell.Np, obj.mesh.K, obj.NonhydroRHS );
            
            Index = physClass.meshUnion.z == 0;
            NonhydroPressure(Index) = 0;
            
%             Index = abs(NonhydroPressure)<1e-10;
%             NonhydroPressure(Index) = 0;
            
            
            %checked
            fphys{1}(:,:,obj.varIndex(1:3)) = mxUpdateConservativeFinalVelocity( NonhydroPressure, fphys{1}, obj.varIndex, ...
                obj.rho, deltatime, obj.PSPX, obj.PSPY, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge,...
                obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype), obj.mesh.EToE);
            
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
        
        function matGetBottomOldMomentum( obj, physClass, fphys )
            edge = physClass.meshUnion.BottomBoundaryEdge;
            [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            obj.BotBEUold = fm(:,:,obj.varIndex(1));
            obj.BotBEVold = fm(:,:,obj.varIndex(2));
            obj.BotBEWold = fm(:,:,obj.varIndex(3));
            edge = physClass.meshUnion.BoundaryEdge;
            [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            obj.BEUold = fm(:,:,obj.varIndex(1));
            obj.BEVold = fm(:,:,obj.varIndex(2));
            obj.BEWold = fm(:,:,obj.varIndex(3));
        end
        
        function matGetBottomNewMomentum( obj, physClass, fphys )
            edge = physClass.meshUnion.BottomBoundaryEdge;
            [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            obj.BotBEUnew = fm(:,:,obj.varIndex(1));
            obj.BotBEVnew = fm(:,:,obj.varIndex(2));
            obj.BotBEWnew = fm(:,:,obj.varIndex(3));
            edge = physClass.meshUnion.BoundaryEdge;
            [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            obj.BEUnew = fm(:,:,obj.varIndex(1));
            obj.BEVnew = fm(:,:,obj.varIndex(2));
            obj.BEWnew = fm(:,:,obj.varIndex(3));
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
            
            %checked
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, ~, ~, ~, obj.PHPX, obj.PHPY ] = ...
                mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            
            %              Rhs = -fphys{1}(:,:,obj.varIndex(4)) .* ( obj.PUPX + obj.PVPY + obj.PUPS .* obj.PSPX + obj.PVPS .* obj.PSPY );
            %
            %              fphys{1}(:,:,obj.varIndex(3)) = fphys{1}(:,:,obj.varIndex(4)) .* physClass.meshUnion(1).VerticalIntegralField(Rhs);
            
            %             fphys{1}(:,:,obj.varIndex(3)) = mxUpdateVerticalFinalVelocity( fphys{1}, obj.varIndex, physClass.hcrit, obj.PSPX,...
            %                 obj.PSPY, obj.PUPX, obj.PVPY, obj.PUPS, obj.PVPS, physClass.VerticalVelocitySolver.RHSCoeMatrix{1},...
            %                 physClass.VerticalVelocitySolver.VertCoeMatrix{1}, obj.mesh, obj.cell, obj.BottomBoundaryEdge);
            
            fphys{1}(:,:,obj.varIndex(3)) = mxUpdateVerticalFinalVelocityIntegralForm( fphys{1}, obj.varIndex, physClass.hcrit, obj.PSPX,...
                obj.PSPY, obj.PUPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.mesh, obj.cell, obj.BottomBoundaryEdge);
        end
        
        function matClearGlobalMemory( obj )
            clear mxCalculatePartialDerivativeUpdated;
            clear mxAssembleGlobalStiffMatrixNew;
            clear mxImposeBoundaryCondition;
            clear mxCalculateBottomVerticalVelocity;
            clear mxNonhydroSystemSolve;
        end
        
        % checked
        function TestPartialDerivativeCalculation(obj, physClass, fphys, fphys2d)
            
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew, obj.Unew, obj.Vnew, obj.PHPX, obj.PHPY ] = mxCalculatePartialDerivativeUpdated( physClass.hcrit,...
                obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            
        end
        
        %checked
        function RHS = TestAssembleRHS( obj, physClass, fphys, deltatime, PUPX, PUPS, PVPY, PVPS, PWPS, PSPX, PSPY )
            
            RHS = mxAssembleNonhydroRHS(PUPX, PUPS, PVPY, PVPS, PWPS, PSPX, ...
                PSPY, fphys{1}(:,:,obj.varIndex(4)), deltatime, obj.rho, physClass.hcrit, obj.mesh.J, obj.cell.M);
            
        end
        
        function Matrix = TestAssembleGlobalStiffMatrix( obj, physClass, fphys, fphys2d )
            %checked
            [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew, obj.Unew, obj.Vnew, obj.PHPX, obj.PHPY ] = ...
                mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            
            Matrix = mxAssembleGlobalStiffMatrixNew(obj.PSPX, obj.PSPY, obj.SQPSPX, ...
                obj.SQPSPY, physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.InnerEdge, obj.BottomEdge,...
                obj.SurfaceBoundaryEdge, obj.cell2d.M, obj.mesh2d.J, obj.mesh2d.K, physClass.SurfaceBoundaryEdgeType,...
                obj.SortedIEnx, obj.SortedIEny, obj.SortedIEGlobalFace, obj.SortedIEAdjEle, obj.SortedIEReverseFlag,...
                obj.SortedIEInternalFace);
            
        end
        
        function fphys = TestUpdateConservativeFinalVelocity( obj, physClass, NonhydroPressure, fphys, deltatime, PSPX, PSPY )
            
            fphys{1}(:,:,obj.varIndex(1:3)) = mxUpdateConservativeFinalVelocity( NonhydroPressure, fphys{1}, obj.varIndex, ...
                obj.rho, deltatime, PSPX, PSPY, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge,...
                obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype));
        end
        
        %Checked
        function TestNewFormGlobalStiffMatrix( obj, physClass, fphys, fphys2d )
            % The coefficient is constant
            if strcmp(physClass.ConstantCoe, 'True')
                obj.PSPX = physClass.K13;
                obj.PSPY = physClass.K23;
                obj.SQPSPX = obj.PSPX .* obj.PSPX;
                obj.SQPSPY = obj.PSPY .* obj.PSPY;
            else
                [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                    obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew, obj.Unew, obj.Vnew, obj.PHPX, obj.PHPY ] = ...
                    mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                    obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                    int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                    fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                    int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            end
            
            obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrixNew(obj.SPNPX, obj.SPNPY, obj.PSPX, obj.PSPY, obj.SQPSPX, ...
                obj.SQPSPY, physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.InnerEdge, obj.BottomEdge,...
                obj.SurfaceBoundaryEdge, obj.cell2d.M, ...
                obj.mesh2d.J, obj.mesh2d.K, physClass.SurfaceBoundaryEdgeType,...
                obj.SortedIEnx, obj.SortedIEny, obj.SortedIEGlobalFace, obj.SortedIEAdjEle, obj.SortedIEReverseFlag,...
                obj.SortedIEInternalFace);
            
        end
        
        %checked
        function [ StiffMatrix, RHS] = TestImposeBoundaryCondition( obj, physClass, fphys, StiffMatrix, RHS, PSPX, PSPY, Wold, Wnew, deltatime, ...
                PWPS, Unew, Uold, Vnew, Vold, PUPX, PUPY, PUPS, PVPX, PVPY, PVPS, PHPX, PHPY)
            
            [ StiffMatrix, RHS] = mxImposeBoundaryCondition( StiffMatrix, PSPX, PSPY, ...
                physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.BottomBoundaryEdge,...
                obj.BoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype), obj.mesh2d, obj.cell2d, obj.InnerEdge2d,...
                obj.BoundaryEdge2d, Wold, Wnew, deltatime, obj.rho, fphys{1}(:,:,obj.varIndex(1)), ...
                fphys{1}(:,:,obj.varIndex(2)), RHS, PWPS, obj.BoundNonhydroPressure, Unew, Uold, Vnew, Vold,...
                PUPX, PUPY, PUPS, PVPX, PVPY, PVPS, PHPX, PHPY, physClass.gra, fphys{1}(:,:,obj.varIndex(5)), ...
                fphys{1}(:,:,obj.varIndex(6)), obj.BoundNonhydroGrad);
            
        end
        
        % Checked
        function TestAssembleFinalGlobalStiffMatrix( obj, physClass, fphys, StiffMatrix, PHPX, PHPY )
            
            obj.GlobalStiffMatrix = mxAssembleFinalGlobalStiffMatrix(obj.cell.Np, obj.mesh.K, physClass.hcrit, obj.mesh.EToE, obj.cell.Nface,...
                obj.cell.M, obj.mesh.J, StiffMatrix, obj.PNPX, obj.PNPY, obj.PNPS, fphys{1}(:,:,obj.varIndex(4)), ...
                PHPX, PHPY, fphys{1}(:,:,obj.varIndex(5)), fphys{1}(:,:,obj.varIndex(6)), obj.mesh.z, obj.UniEleNumber, obj.UniEle);
        end
        
        function AssembleStiffMatrix( obj, physClass, fphys, fphys2d )
            
            % The coefficient is constant
            if strcmp(physClass.ConstantCoe, 'True')
                obj.PSPX = physClass.K13;
                obj.PSPY = physClass.K23;
                obj.SQPSPX = obj.PSPX .* obj.PSPX;
                obj.SQPSPY = obj.PSPY .* obj.PSPY;
                [ ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj.PHPX, obj.PHPY ] = ...
                    mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                    obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                    int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                    fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                    int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            else
                [ obj.PSPX, obj.PSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, obj.PUPY, ...
                    obj.PVPX, obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS, obj.Wnew, obj.Unew, obj.Vnew, obj.PHPX, obj.PHPY ] = ...
                    mxCalculatePartialDerivativeUpdated( physClass.hcrit, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, ...
                    obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, ...
                    int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),...
                    fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                    int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), physClass.frhs2d{1}(:,:,1));
            end
            
            obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrixNew(obj.SPNPX, obj.SPNPY, obj.PSPX, obj.PSPY, obj.SQPSPX, ...
                obj.SQPSPY, physClass.hcrit, fphys{1}(:,:,obj.varIndex(4)), obj.mesh, obj.cell, obj.InnerEdge, obj.BottomEdge,...
                obj.SurfaceBoundaryEdge, obj.cell2d.M, ...
                obj.mesh2d.J, obj.mesh2d.K, physClass.SurfaceBoundaryEdgeType,...
                obj.SortedIEnx, obj.SortedIEny, obj.SortedIEGlobalFace, obj.SortedIEAdjEle, obj.SortedIEReverseFlag,...
                obj.SortedIEInternalFace);
            
            obj.GlobalStiffMatrix = mxAssembleFinalGlobalStiffMatrix(obj.cell.Np, obj.mesh.K, physClass.hcrit, obj.mesh.EToE, obj.cell.Nface,...
                obj.cell.M, obj.mesh.J, obj.GlobalStiffMatrix, obj.PNPX, obj.PNPY, obj.PNPS, fphys{1}(:,:,obj.varIndex(4)), ...
                obj.PHPX, obj.PHPY, fphys{1}(:,:,obj.varIndex(5)), fphys{1}(:,:,obj.varIndex(6)), obj.mesh.z, obj.UniEleNumber, obj.UniEle);
        end
        
    end
    
end

