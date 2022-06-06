classdef NdgQuadFreeStrongFormAdvSWE3dSolver3d < NdgQuadFreeStrongFormAdvSolver3d
    
    properties
        mesh
        cell
        InnerEdge
        BoundaryEdge
        BottomEdge
        BottomBoundaryEdge
        SurfaceBoundaryEdge
    end
    
    methods
        function obj = NdgQuadFreeStrongFormAdvSWE3dSolver3d( phys )
            obj = obj@NdgQuadFreeStrongFormAdvSolver3d( phys );
            warning('off');
            obj.mesh = struct(phys.meshUnion);
            obj.cell = struct(phys.meshUnion.cell);
            obj.InnerEdge = struct(phys.meshUnion.InnerEdge);
            obj.BoundaryEdge = struct(phys.meshUnion.BoundaryEdge);
            obj.BottomEdge = struct(phys.meshUnion.BottomEdge);
            obj.BottomBoundaryEdge = struct(phys.meshUnion.BottomBoundaryEdge);
            obj.SurfaceBoundaryEdge = struct(phys.meshUnion.SurfaceBoundaryEdge);
            warning('on');
            obj.matClearGlobalMemory( );
        end
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, physClass, fphys )
            % evaluate inner edge
            % %             tic;
            %            [ physClass.frhs{1}, TempE, TempG, TempH, TempVolumeIntegral,...
            %                OutputIEFluxM, OutputIEFluxP, OutputIEFluxS, OutputBEFluxM, OutputBEFluxP, OutputBEFluxS, OutputBotEFluxM, ...
            %                OutputBotEFluxP, OutputBotEFluxS, OutputBotBEFluxM, OutputBotBEFluxS, OutputSurfBEFluxM,...
            %                OutputSurfBEFluxS ]= mxEvaluateQuadFreeStrongFormAdvSWE3dRHS( obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge,...
            %                 obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, physClass.varFieldIndex, fphys{1}, ...
            %                  physClass.fext3d{1}, physClass.gra, physClass.hcrit, int8(physClass.meshUnion.BoundaryEdge.ftype), physClass.SurfaceDate);
            
            [ physClass.frhs{1} ]= mxEvaluateQuadFreeStrongFormAdvSWE3dRHS( obj.mesh, obj.cell, obj.InnerEdge, obj.BoundaryEdge,...
                obj.BottomEdge, obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, physClass.varFieldIndex, fphys{1}, ...
                physClass.fext3d{1}, physClass.gra, physClass.hcrit, int8(physClass.meshUnion.BoundaryEdge.ftype) );
            
            
%             disp('The maximum difference for the rhs3d between Matlab version and C version is:')
%             disp(max(max(-(physClass.meshUnion.rx .* ( physClass.meshUnion.cell.Dr * (fphys{1}(:, :, 1).*fphys{1}(:, :, 14)./fphys{1}(:, :, 4)) ) + physClass.meshUnion.sx .* ( physClass.meshUnion.cell.Ds * (fphys{1}(:, :, 1).*fphys{1}(:, :, 14)./fphys{1}(:, :, 4)) ) +...
%                 physClass.meshUnion.ry .* ( physClass.meshUnion.cell.Dr * (fphys{1}(:, :, 2).*fphys{1}(:, :, 14)./fphys{1}(:, :, 4)) ) + physClass.meshUnion.sy .* ( physClass.meshUnion.cell.Ds * (fphys{1}(:, :, 2).*fphys{1}(:, :, 14)./fphys{1}(:, :, 4)) )) - physClass.frhs{1}(:,:,3))));            
            
            % %             Cversion = toc;
            % %             tic;
            % %             for m = 1:physClass.Nmesh
            % %                 mesh3d = physClass.meshUnion(m);
            % %
            % %                 edge = mesh3d.InnerEdge;
            % %                 [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
            % %
            % %                 [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
            % %                 IEFluxM = fluxM(:,:,[2,3]);
            % %                 [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
            % %                 IEFluxP = fluxP(:,:,[2,3]);
            % %                 [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), fp(:,:,[4, 1, 2]), edge );
            % %                 IEFluxS = fluxS(:,:,[2,3]);
            % %                 [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM(:,:,[2,3]), fluxP(:,:,[2,3]), fluxS(:,:,[2,3]) );
            % %
            % %                 edge = mesh3d.BoundaryEdge;
            % %                 [ fm, physClass.BoundaryEdgefp{m} ] = edge.matEvaluateSurfValue( fphys );
            % %
            % %                 [ fm, physClass.BoundaryEdgefp{m} ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, fm, physClass.BoundaryEdgefp{m}, physClass.fext3d{m} );
            % %                 [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
            % %                 BEFluxM = fluxM(:,:,[2,3]);
            % %                 [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), physClass.BoundaryEdgefp{m}(:,:,[4, 1, 2]), edge );
            % %                 BEFluxS = fluxS(:,:,[2,3]);
            % %                 [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), fluxS(:,:,[2,3]) );
            % %
            % %
            % %                 edge = mesh3d.BottomEdge;
            % %                 [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
            % %                 [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
            % %                 BotEFluxM = fluxM(:,:,[2,3]);
            % %                 [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
            % %                 BotEFluxP = fluxP(:,:,[2,3]);
            % %                 [ OmegafluxS(:,:,1) ] = 0.5*edge.nz.*( fm(:,:,1).*fm(:,:,3)./fm(:,:,4) + fp(:,:,1).*fp(:,:,3)./fp(:,:,4) );
            % %                 [ OmegafluxS(:,:,2) ] = 0.5*edge.nz.*( fm(:,:,2).*fm(:,:,3)./fm(:,:,4) + fp(:,:,2).*fp(:,:,3)./fp(:,:,4) );
            % %                 BotEFluxS = OmegafluxS;
            % % %                 [ OmegafluxS(:,:,1) ] = edge.nz.*( ( fp(:,:,3)>=0 ) .* fp(:,:,1).*fp(:,:,3)./fp(:,:,4) + ( fp(:,:,3) < 0 ) .* fm(:,:,1).*fp(:,:,3)./fm(:,:,4) ) ;
            % % %                 [ OmegafluxS(:,:,2) ] = edge.nz.*( ( fp(:,:,3)>=0 ) .* fp(:,:,2).*fp(:,:,3)./fp(:,:,4) + ( fp(:,:,3) < 0 ) .* fm(:,:,2).*fp(:,:,3)./fm(:,:,4) ) ;
            % %                 [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), fluxP(:,:,[2,3]), OmegafluxS );
            % %
            % %                 edge = mesh3d.BottomBoundaryEdge;
            % %                 [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            % %                 [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
            % %                 BotBEFluxM = fluxM(:,:,[2,3]);
            % %                 [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), zeros(size(fluxM(:,:,[2,3]) )) );
            % %                 BotBEFluxS = zeros(size(fluxM(:,:,[2,3]) ));
            % %
            % %                 edge = mesh3d.SurfaceBoundaryEdge;
            % %                 [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
            % %                 [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
            % %                 SurfBEFluxM = fluxM(:,:,[2,3]);
            % %                 [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), zeros(size(fluxM(:,:,[2,3]) )) );
            % %                 SurfBEFluxS = zeros(size(fluxM(:,:,[2,3]) ));
            % %             end
            % %
%             for m = 1:physClass.Nmesh % calculate RHS term on each mesh
%                 mesh = physClass.meshUnion(m);
%                 [ E, G, H ] = physClass.matEvaluateFlux( mesh, fphys{m} );
%                 Tempfrhs = zeros(mesh.cell.Np, mesh.K, 4);
%                 for i = 1:2
%                     %                     physClass.frhs{m}(:,:,i) = ...
%                     %                         physClass.frhs{m}(:,:,i) ...
%                     %                         - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
%                     %                         - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
%                     %                         - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
%                     %                         - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ) ...
%                     %                         - obj.tz{m}.*( obj.Dt{m} * H(:,:,i) );
%                     Tempfrhs(:,:,i) = ...
%                         - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
%                         - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
%                         - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
%                         - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ) ...
%                         - obj.tz{m}.*( obj.Dt{m} * H(:,:,i) );
%                 end
%                 E(:,:,1) = fphys{m}(:,:,1) .* fphys{m}(:,:,14)./fphys{m}(:,:,4);
%                 G(:,:,1) = fphys{m}(:,:,2) .* fphys{m}(:,:,14)./fphys{m}(:,:,4);
%                 E(:,:,1) = fphys{m}(:,:,1) .* 10;
%                 G(:,:,1) = fphys{m}(:,:,2) .* 10;
%                 H(:,:,1) = zeros(size(E(:,:,1)));
%                 E(:,:,2) = fphys{m}(:,:,1) .* fphys{m}(:,:,15)./fphys{m}(:,:,4);
%                 G(:,:,2) = fphys{m}(:,:,2) .* fphys{m}(:,:,15)./fphys{m}(:,:,4);
%                 H(:,:,2) = fphys{m}(:,:,3) .* fphys{m}(:,:,15)./fphys{m}(:,:,4);
% %                 
%                 for i = 3:4
%                     Tempfrhs(:,:,i) = ...
%                         - obj.rx{m}.*( obj.Dr{m} * E(:,:,i-2) ) ...
%                         - obj.sx{m}.*( obj.Ds{m} * E(:,:,i-2) ) ...
%                         - obj.ry{m}.*( obj.Dr{m} * G(:,:,i-2) ) ...
%                         - obj.sy{m}.*( obj.Ds{m} * G(:,:,i-2) ) ...
%                         - obj.tz{m}.*( obj.Dt{m} * H(:,:,i-2) );
%                 end
% %                 
%             end
%             disp("For the three-dimensional part, the maximum difference is:");
%             disp(max(max(abs(physClass.frhs{1}(:,:,3) - Tempfrhs(:,:,3)))));
            % %             matVersion = toc;
            % %             fprintf("The speed up ratio is:%f\n",matVersion/Cversion);
        end
        
        function matClearGlobalMemory( obj )
            clear mxEvaluateQuadFreeStrongFormAdvSWE3dRHS;
        end
    end
    
end