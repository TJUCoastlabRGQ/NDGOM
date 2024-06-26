classdef SWEQuadFreeStrongFormPCESolver2d
    %SWEQUADFREESTRONGFORMPCESOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        mesh
        InnerEdge
        BoundaryEdge
        cell
        mesh3d
        cell3d
        InnerEdge3d
        BoundaryEdge3d        
    end
    
    methods
        
        function obj = SWEQuadFreeStrongFormPCESolver2d( mesh )
            warning('off');
            obj.mesh = struct(mesh.mesh2d);
            obj.InnerEdge = struct(mesh.mesh2d.InnerEdge);
            obj.BoundaryEdge = struct(mesh.mesh2d.BoundaryEdge);
            obj.cell = struct(mesh.mesh2d.cell);
            obj.mesh3d = struct(mesh);
            obj.InnerEdge3d = struct(mesh.InnerEdge);
            obj.BoundaryEdge3d = struct(mesh.BoundaryEdge);
            obj.cell3d = struct(mesh.cell);            
            warning('on');
            obj.matClearGlobalMemory( );
        end
        function evaluateAdvectionRHS( obj, physClass, fphys2d, fphys )
% %             tic;
%             physClass.frhs2d{1} = mxEvaluatePCERHS( obj.mesh, obj.InnerEdge, obj.BoundaryEdge,...
%                 obj.cell, int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), fphys2d{1}, physClass.fext2d{1}, physClass.gra, physClass.hcrit );
            
            physClass.frhs2d{1} = mxEvaluatePCERHSUpdated( obj.mesh, obj.mesh3d, obj.cell3d, obj.InnerEdge, obj.BoundaryEdge,...
                obj.cell, int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), int8(physClass.meshUnion.BoundaryEdge.ftype),...
                fphys2d{1}, fphys{1}, physClass.fext2d{1}, physClass.fext3d{1}, physClass.gra, physClass.hcrit, ...
                obj.InnerEdge3d, obj.BoundaryEdge3d);
            
%             disp('The maximum difference for the rhs2d and the integral of rhs3d is:');
% %             disp(max(max(abs(physClass.frhs2d{1} - physClass.meshUnion.VerticalColumnIntegralField(physClass.frhs{1}(:,:,3))./10))));
%             disp(max(max(abs((physClass.frhs2d{1} - physClass.meshUnion.VerticalColumnIntegralField(physClass.frhs{1}(:,:,3))./10)./physClass.frhs2d{1}) * 100)));
%             
%             disp('The maximum difference for the rhs2d between Matlab version and C version is:')
%             disp(max(max(-(physClass.mesh2d.rx .* ( physClass.mesh2d.cell.Dr * fphys2d{1}(:, :, 2) ) + physClass.mesh2d.sx .* ( physClass.mesh2d.cell.Ds * fphys2d{1}(:, :, 2) ) +...
%                 physClass.mesh2d.ry .* ( physClass.mesh2d.cell.Dr * fphys2d{1}(:, :, 3) ) + physClass.mesh2d.sy .* ( physClass.mesh2d.cell.Ds * fphys2d{1}(:, :, 3) )) - physClass.frhs2d{1})));
            
% %             t1 = toc;
% %             tic;
            % evaluate inner edge
%             for m = 1:physClass.Nmesh
%                 
%                 %  Function used to calculate the vertically averaged horizontal momentum term
%                 %                 mesh3d = physClass.meshUnion(m);
%                 mesh2d = physClass.mesh2d(m);
%                 
%                 %Function used to calculate the two dimentional PCE volume term
%                 Tempfrhs2d = -( ...
%                     mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d{m}(:, :, 2) ) + ...
%                     mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d{m}(:, :, 2) ) + ...
%                     mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d{m}(:, :, 3) ) + ...
%                     mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d{m}(:, :, 3) ) );
% 
%                 Tempfrhs2d = -( ...
%                     physClass.meshUnion.VerticalColumnIntegralField(physClass.meshUnion.rx .* ( physClass.meshUnion.cell.Dr * fphys{m}(:, :, 1) ) + ...
%                     physClass.meshUnion.sx .* ( physClass.meshUnion.cell.Ds * fphys{m}(:, :, 1) ) + ...
%                     physClass.meshUnion.ry .* ( physClass.meshUnion.cell.Dr * fphys{m}(:, :, 2) ) + ...
%                     physClass.meshUnion.sy .* ( physClass.meshUnion.cell.Ds * fphys{m}(:, :, 2) ) ));

%                 Edge = mesh2d.InnerEdge;
%                 [fm, fp] = Edge.matEvaluateSurfValue( fphys2d );
%                 [ fluxM ] = obj.matEvaluateSurfFlux( Edge, Edge.nx, Edge.ny, fm );
%                 [ fluxP ] = obj.matEvaluateSurfFlux( Edge, Edge.nx, Edge.ny, fp );
%                 [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh2d, Edge.nx, Edge.ny, fm, fp, Edge);
%                 physClass.frhs2d{m} = physClass.frhs2d{m} + Edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );
%                 
%                 Edge = mesh2d.BoundaryEdge;
%                 [fm, fp] = Edge.matEvaluateSurfValue( fphys2d );
%                 [fm, fp] = obj.matImposeBoundaryCondition( Edge, Edge.nx, Edge.ny, fm, fp, physClass.fext2d{m});
%                 [ fluxM ] = obj.matEvaluateSurfFlux( Edge, Edge.nx, Edge.ny, fm );
%                 [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh2d, Edge.nx, Edge.ny, fm, fp, Edge);
%                 physClass.frhs2d{m} = physClass.frhs2d{m} + Edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                
                % Function used to calculate the two dimentional PCE inner surface term
                %                 InnerEdge = mesh3d.InnerEdge;
                %                 InnerEdge2d = mesh2d.InnerEdge;
                %                 physClass.InnerEdgeFluxM2d{m} = InnerEdge.VerticalColumnIntegralField( physClass.InnerEdgeFluxM{m}(:,:,1) );
                %                 physClass.InnerEdgeFluxP2d{m} = InnerEdge.VerticalColumnIntegralField( physClass.InnerEdgeFluxP{m}(:,:,1) );
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
                %                 physClass.InnerEdgeFluxS2d{m} = InnerEdge.VerticalColumnIntegralField( physClass.InnerEdgeFluxS{m}(:,:,1) );
                %                 physClass.frhs2d{m} = physClass.frhs2d{m} + InnerEdge2d.matEvaluateStrongFormEdgeRHS( physClass.InnerEdgeFluxM2d{m}, physClass.InnerEdgeFluxP2d{m}, physClass.InnerEdgeFluxS2d{m} );
                
                % Function used to calculate the two dimentional PCE boundary surface integration term
                %                 BoundaryEdge = mesh3d.BoundaryEdge;
                %                 BoundaryEdge2d = mesh2d.BoundaryEdge;
                
                %                 physClass.BoundaryEdgeFluxM2d{m} = BoundaryEdge.VerticalColumnIntegralField( physClass.BoundaryEdgeFluxM{m}(:,:,1) );
                %                 physClass.BoundaryEdgeFluxS2d{m} = BoundaryEdge.VerticalColumnIntegralField( physClass.BoundaryEdgeFluxS{m}(:,:,1) );
                
                %                 physClass.frhs2d{m} = physClass.frhs2d{m} + BoundaryEdge2d.matEvaluateStrongFormEdgeRHS( physClass.BoundaryEdgeFluxM2d{m}, physClass.BoundaryEdgeFluxS2d{m} );
                
%             end
% %             t2 = toc;
% %             fprintf("The speed up ratio is:%f\n",t2/t1);
%            disp("For the rhs2d, the maxmum difference is:");
%            disp(max(max(abs((10*Tempfrhs2d - physClass.meshUnion.VerticalColumnIntegralField(physClass.frhs{1}(:,:,3)))))));
        end
        
        function matClearGlobalMemory( obj )
            clear mxEvaluatePCERHS;
            clear mxEvaluatePCERHSUpdated;
        end
        
    end
end

