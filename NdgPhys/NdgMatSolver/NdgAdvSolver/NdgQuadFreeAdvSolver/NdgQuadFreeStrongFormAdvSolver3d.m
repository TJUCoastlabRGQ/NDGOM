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
                mesh3d = phys.mesh3d(m);
                edge = mesh3d.InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, edge.ny, edge.nz, fm, fp );
                [ phys.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );

                edge = mesh3d.BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = phys.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, fm, fp, phys.fext3d{m} );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, edge.ny, edge.nz, fm, fp );
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                edge = mesh3d.SurfaceBoundaryEdge;
                [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                % At present, this solver is designed for three dimensional
                % barotropic shallow water equation, and the top boundary
                % condition is not considered here, and I simply set this
                % numerical flux to be zero, this would be corrected later
                [ fluxS ] =  zeros(size(fluxM));
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                edge = mesh3d.BottomEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, edge.ny, edge.nz, fm, fp );
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );
                
                edge = mesh3d.BottomBoundaryEdge;
                [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                % At present, this solver is designed for three dimensional
                % barotropic shallow water equation, and the top boundary
                % condition is not considered here, and I simply set this
                % numerical flux to be zero, this would be corrected later                
                [ fluxS ] =  zeros(size(fluxM));
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );    
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