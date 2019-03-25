function matEvaluateRHS( obj, fphys2d, fphys3d )
%MATEVALUATERHS Summary of this function goes here
%   Detailed explanation goes here

for m = 1:obj.Nmesh
    mesh2d = obj.mesh2d(m);
    mesh3d = obj.mesh3d(m);
    
    % evaluate 2d PCE volume integral term
    fphys2d{m} = obj.matEvaluate2dHorizonMomentum( mesh3d, ...
        fphys2d{m}, fphys3d{m} );
    
    obj.frhs2d{m} = obj.matEvaluate2dHorizonPCEVolumeTerm( mesh2d, fphys2d{m} );

    % evaluate 2d PCE surface integral term
    obj.frhs2d{m} = obj.frhs2d{m} + obj.matEvaluate2dHorizonPCEInnerSurfaceTerm( ...
        mesh2d.InnerEdge, fphys2d );
    
    obj.frhs2d{m} = obj.frhs2d{m} + obj.matEvaluate2dHorizonPCEBoundaryTerm( ...
         mesh2d.BoundaryEdge, fphys2d, obj.fext2d{m});
    
    % extend 3d field
    fphys3d = obj.matEvaluate3dAuxiliaryVariable(  mesh3d, fphys2d, fphys3d );   %1
    
    % evaluate 3d velocity volume integral
    obj.frhs3d{m} = obj.matEvaluate3dVolumeTerm( mesh3d, fphys3d{m} );   %2
    
    % evaluate 3d velocity horizontal surface integral
    obj.frhs3d{m} = obj.frhs3d{m} +obj.matEvaluate3dSideSurfaceTerm( ...           %3
         mesh3d.InnerEdge, fphys3d );
    
    obj.frhs3d{m} = obj.frhs3d{m} + obj.matEvaluate3dHorizontalBoundaryTerm( ...     %4
         mesh3d.BoundaryEdge, fphys3d, obj.fext3d{m} );
    
    % evaluate 3d velocity field surface integral
    obj.frhs3d{m} = obj.frhs3d{m} + obj.matEvaluate3dSurfaceBoundaryTerm( ...        %5
         mesh3d.SurfaceBoundaryEdge, fphys3d );
    
    % evaluate 3d velocity field bottom integral
    obj.frhs3d{m} = obj.frhs3d{m} + obj.matEvaluate3dBottomEdgeTerm( ...         %6
         mesh3d.BottomEdge, fphys3d);
    
    obj.frhs3d{m} = obj.frhs3d{m} + obj.matEvaluate3dBottomBoundaryTerm( ...       %7
         mesh3d.BottomBoundaryEdge, fphys3d );
    
    obj.frhs3d{m} = obj.frhs3d{m} + obj.matEvaluate3dVerticalAuxialaryVariableFaceTerm( ...       %7
         mesh3d, fphys3d );    
end

end
