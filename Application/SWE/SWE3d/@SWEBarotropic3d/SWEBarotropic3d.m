classdef SWEBarotropic3d < SWEAbstract3d
    %SWEBAROTROPIC3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    properties ( Constant )
        %> number of physical field
        Nfield2d = 7 %[eta HU HV H Z Zx Zy]
        %> num of 3d physical field
        Nfield3d = 11 % [ Hu, Hv, w, tau_x, tau_y, H, eta, Z, Zx, Zy, MUy ]
        %> number of variable field
        Nvar = 3
        %> index of variable in physical field
        varFieldIndex = [ 1, 2 ]
    end
    
    
    methods
        
        function matSolve( obj )
            matEvaluateRK45( obj );
        end
        
        matEvaluateRK45( obj );
    end
    
    methods ( Access = protected )
        
        %         Volume_frhs3d = matEvaluate3dVolumeTerm( obj, mesh3d, fphys3d );
        %
        %         InnerSurface_frhs3d = matEvaluate3dInnerSurfaceTerm( obj, InnerEdge, fphys3d );
        %
        %         VerticalBoundarySurface_frhs3d = matEvaluate3dVerticalBoundarySurfaceTerm( obj, BoundaryEdge, fphys3d );
        %
        %         Bottom_frhs3d = matEvaluate3dBottomTerm(obj, BottomEdge, fphys3d);
        %
        %         SurfaceBoundary_frhs3d =  matEvaluate3dSurfaceBoundaryTerm( obj, SurfaceBoundaryEdge, fphys3d );
        %
        %         BottomBoundary_frhs3d =  matEvaluate3dBottomBoundaryTerm( obj, BottomBoundaryEdge, fphys3d );
        
    end
    
end

