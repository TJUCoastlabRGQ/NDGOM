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
        
        NonhydroRHS
        
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
    
    methods
        function obj = NdgQuadratureFreeNonhydrostaticSolver3d( PhysClass, mesh )
            obj.matSetInitializeCharacteristicMatrix( mesh );
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
            
            [ obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS ] = mxCalculatePartialDerivative( physClass.hcrit,...
                  obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                  obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                  physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                 obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                 int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));      
             
             obj.GlobalStiffMatrix = mxAssembleGlobalStiffMatrix( obj.PNPS, obj.SPNPX, obj.SPNPY, obj.SPNPS, ...
                 obj.MSPNPX, obj.MSPNPY, obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY,...
                 fphys{1}(:,:,obj.varIndex(4)), physClass.hcrit, obj.BoundaryEdge2d, ...
                 obj.mesh, obj.cell, obj.BoundaryEdge, int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));
             
             obj.NonhydroRHS = mxAssembleNonhydroRHS(obj.PUPX, obj.PUPS, obj.PVPY, obj.PVPS, obj.PWPS, obj.PSPX, ...
                 obj.PSPY, fphys{1}(:,:,obj.varIndex(4)), deltatime, obj.rho, physClass.hcrit);
             
             NonhydroPressure = obj.GlobalStiffMatrix\obj.NonhydroRHS;
                          
             fphys{1}(:,:,obj.varIndex(1:3)) = mxUpdateConservativeFinalVelocity( NonhydroPressure, fphys{1}, obj.varIndex, ...
                 obj.rho, deltatime, obj.PSPX, obj.PSPY, obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge,...
                 obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, int8(physClass.meshUnion.BoundaryEdge.ftype));
        end  
        
        function matClearGlobalMemory( obj )
            clear mxCalculatePartialDerivative;
        end
        
        function TestPartialDerivativeCalculation(obj, physClass, fphys, fphys2d, ~)
            
             [ obj.PSPX, obj.PSPY, obj.SPSPX, obj.SPSPY, obj.SQPSPX, obj.SQPSPY, obj.PUPX, ...
                obj.PVPY, obj.PUPS, obj.PVPS, obj.PWPS ] = mxCalculatePartialDerivative( physClass.hcrit,...
                  obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge, obj.BottomEdge, obj.BottomBoundaryEdge, ...
                  obj.SurfaceBoundaryEdge, fphys{1}, obj.varIndex, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
                  physClass.gra, physClass.fext3d{ 1 }, fphys2d{1}(:,:,1),  fphys2d{1}(:,:,4), physClass.fext2d{ 1 }, ...
                 obj.mesh2d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.cell2d, ...
                 int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype));  
             
        end
        
        
    end
    
end

