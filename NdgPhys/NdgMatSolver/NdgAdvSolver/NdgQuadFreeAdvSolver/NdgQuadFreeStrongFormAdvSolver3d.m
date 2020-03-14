classdef NdgQuadFreeStrongFormAdvSolver3d < NdgQuadFreeStrongFormSolver & ...
        NdgAbstractAdvSolver
    
    methods
        function obj = NdgQuadFreeStrongFormAdvSolver3d( phys )
            obj = obj@NdgQuadFreeStrongFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver(phys);
        end
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            
            % evaluate inner edge
            for m = 1:phys.Nmesh
                mesh3d = phys.meshUnion(m);
                edge = mesh3d.InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), fp(:,:,[4, 1, 2]), edge );
                [ phys.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS(:,:,[2,3]) );

                edge = mesh3d.BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = phys.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, fm, fp, phys.fext3d{m} );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), fp(:,:,[4, 1, 2]), edge );
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS(:,:,[2,3]) );
                
                % we note that for the three dimensional nonlinear shallow
                % water equation, Newmann boundary about the velocity is
                % imposed at both the surface and bottom boundary, under
                % this kind of circumastances, the numerical flux is the
                % same with the local numerical flux, so the contribution
                % of these term to the right hand is canceled for the
                % advective part
                
%                 edge = mesh3d.SurfaceBoundaryEdge;
%                 [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
%                 [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
%                 [ fluxS ] =  zeros(size(fluxM));
%                 [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
%                 
                edge = mesh3d.BottomEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
                [ OmegafluxS(:,:,1) ] = 0.5*edge.nz.*(fm(:,:,1).*fm(:,:,3)./fm(:,:,4)+fp(:,:,1).*fp(:,:,3)./fp(:,:,4));
                [ OmegafluxS(:,:,2) ] = 0.5*edge.nz.*(fm(:,:,2).*fm(:,:,3)./fm(:,:,4)+fp(:,:,2).*fp(:,:,3)./fp(:,:,4));
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, OmegafluxS );
                
%                 edge = mesh3d.BottomBoundaryEdge;
%                 [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
%                 [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );            
% %                 [ fluxS ] =  zeros(size(fluxM));
%                 [ fluxS ] = fluxM;
%                 [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );    
            end
            
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E, G, H ] = phys.matEvaluateFlux( mesh, fphys{m} );
                
                for i = 1:phys.Nvar
                    phys.frhs{m}(:,:,i) = ...
                        phys.frhs{m}(:,:,i) + ...
                        - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
                        - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
                        - obj.tx{m}.*( obj.Dt{m} * E(:,:,i) ) ...
                        - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
                        - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ) ...
                        - obj.ty{m}.*( obj.Dt{m} * G(:,:,i) ) ...
                        - obj.rz{m}.*( obj.Dr{m} * H(:,:,i) ) ...
                        - obj.sz{m}.*( obj.Ds{m} * H(:,:,i) ) ...
                        - obj.tz{m}.*( obj.Dt{m} * H(:,:,i) );
                end
            end
        end
    end
    
end