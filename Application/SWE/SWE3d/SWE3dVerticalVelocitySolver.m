classdef SWE3dVerticalVelocitySolver < handle
    %SWEVERTICALVELOCITYSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        % This corresponds to the matrix applied to the horizontal partial
        % derivative term
        RHSCoeMatrix
        % This corresponds to the matrix applied to the vertical velocity
        % term
        VertCoeMatrix
    end
    
    methods
        function obj = SWE3dVerticalVelocitySolver( mesh2d, mesh3d)
            obj.RHSCoeMatrix = cell(1);
            obj.RHSCoeMatrix{1} = zeros( mesh3d.cell.Np, mesh3d.cell.Np, mesh3d.K );
            obj.VertCoeMatrix = obj.RHSCoeMatrix;
            BotEidM = mesh3d.cell.Fmask( mesh3d.cell.Fmask( :,end-1) ~= 0, end-1 );
            Nz = mesh3d.Nz;
            for i = 1:mesh2d.K
                M3d = zeros(mesh3d.cell.Np);
                M3d( BotEidM, BotEidM ) = diag(mesh2d.J(:,i)) * mesh2d.cell.M;     
                for j = 1 : Nz-1
                    CoeMatrix = diag(mesh3d.J(:, (i-1)*Nz + j)) * mesh3d.cell.M * diag(mesh3d.tz(:, (i-1)*Nz + j)) * mesh3d.cell.Dt + M3d;
                    obj.RHSCoeMatrix{1}(:,:,(i-1)*Nz + j)  = CoeMatrix\(diag(mesh3d.J(:, (i-1)*Nz + j)) * mesh3d.cell.M);
                    obj.VertCoeMatrix{1}(:,:,(i-1)*Nz + j)  = CoeMatrix\M3d;
                end
                    CoeMatrix = diag(mesh3d.J(:, i*Nz)) * mesh3d.cell.M * diag(mesh3d.tz(:, i * Nz)) * mesh3d.cell.Dt + 2 * M3d;
                    obj.RHSCoeMatrix{1}(:,:,(i-1)*Nz + j)  = CoeMatrix\(diag(mesh3d.J(:, (i-1)*Nz + j)) * mesh3d.cell.M);                           
            end
        end
        
        function VerticalVelocity = matCalculateVerticalVelocity( obj, physClass, fphys2d, fphys )
            edge = physClass.meshUnion.InnerEdge;
            edge2d = physClass.mesh2d.InnerEdge;
            [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
            [ FluxM3d ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
            [ FluxM2d ] = edge.VerticalColumnIntegralField(FluxM3d(:,:,1) );
            [ FluxP3d ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
            [ FluxP2d ] = edge.VerticalColumnIntegralField(FluxP3d(:,:,1) );
            [ FluxS3d ] = physClass.matEvaluateSurfNumFlux( physClass.meshUnion, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), fp(:,:,[4, 1, 2]), edge );
            [ FluxS2d ] = edge.VerticalColumnIntegralField(FluxS3d(:,:,1) );
            
            [ InnerSurface_frhs3d ] = edge.matEvaluateStrongFromEdgeRHS( FluxM3d(:,:,1), FluxP3d(:,:,1), FluxS3d(:,:,1) ) ;
            [ InnerSurface_frhs2d ] = edge2d.matEvaluateStrongFormEdgeRHS( FluxM2d(:,:,1), FluxP2d(:,:,1), FluxS2d(:,:,1) ) ;
            
            edge = physClass.meshUnion.BoundaryEdge;
            edge2d = physClass.mesh2d.BoundaryEdge;
            [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
            [ fm, fp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, fm, fp, physClass.fext3d{1} );
            [ FluxM3d ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
            [ FluxM2d ] = edge.VerticalColumnIntegralField( FluxM3d(:,:,1) );
            [ FluxS3d ] = physClass.matEvaluateSurfNumFlux( physClass.meshUnion, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), fp(:,:,[4, 1, 2]), edge );
            [ FluxS2d ] = edge.VerticalColumnIntegralField( FluxS3d(:,:,1) );
            
            [ BoundarySurface_frhs3d ] = edge.matEvaluateStrongFormEdgeRHS( FluxM3d(:,:,1), FluxS3d(:,:,1) ) ;
            [ BoundarySurface_frhs2d ] = edge2d.matEvaluateStrongFormEdgeRHS( FluxM2d(:,:,1), FluxS2d(:,:,1) );
            
            field3d =...
                ( physClass.meshUnion.rx .* (physClass.meshUnion.cell.Dr * fphys{1}(:,:,1)) + physClass.meshUnion.sx .* (physClass.meshUnion.cell.Ds * fphys{1}(:,:,1)) ) - ...
                InnerSurface_frhs3d - BoundarySurface_frhs3d  + ( physClass.meshUnion.ry .* (physClass.meshUnion.cell.Dr * fphys{1}(:,:,2)) + physClass.meshUnion.sy .* (physClass.meshUnion.cell.Ds * fphys{1}(:,:,2)) ) ;
            % Term2d = ;
            field2d = physClass.meshUnion.Extend2dField( physClass.mesh2d.rx .* (physClass.mesh2d.cell.Dr * fphys2d{1}(:,:,2)) + physClass.mesh2d.sx .* (physClass.mesh2d.cell.Ds * fphys2d{1}(:,:,2) ) - ...
                InnerSurface_frhs2d - BoundarySurface_frhs2d + physClass.mesh2d.ry .* (physClass.mesh2d.cell.Dr * fphys2d{1}(:,:,3) ) + physClass.mesh2d.sy .* (physClass.mesh2d.cell.Ds * fphys2d{1}(:,:,3)) );
                        
            Nz = physClass.meshUnion.Nz;
            BotEidM = physClass.meshUnion.cell.Fmask( physClass.meshUnion.cell.Fmask( :,end-1) ~= 0, end-1 );
            UpEidM = physClass.meshUnion.cell.Fmask( physClass.meshUnion.cell.Fmask( :,end) ~= 0, end );
            VerticalVelocity = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
            
            VerticalVelocity(:, Nz:Nz:end) =  permute( sum( bsxfun(@times, obj.RHSCoeMatrix{1}(:,:,Nz:Nz:end), ...
                permute( permute( field2d(:, Nz:Nz:end) - field3d(:, Nz:Nz:end), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
            
            for Layer = 2 : Nz
                BotVertVelocity = zeros( physClass.meshUnion.cell.Np, physClass.mesh2d.K );
                BotVertVelocity( BotEidM, : ) = VerticalVelocity( UpEidM, Nz - Layer + 2 : Nz : end);
                VerticalVelocity(:, Nz - Layer + 1:Nz:end) = permute( sum( bsxfun(@times, obj.RHSCoeMatrix{1}(:,:,Nz - Layer + 1:Nz:end), ...
                permute( permute( field2d(:, Nz - Layer + 1:Nz:end) - field3d(:, Nz - Layer + 1:Nz:end), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]) + permute( sum( bsxfun(@times, obj.VertCoeMatrix{1}(:,:,Nz - Layer + 1:Nz:end), ...
                permute( permute( BotVertVelocity, [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
            end
        end
    end
    
end

