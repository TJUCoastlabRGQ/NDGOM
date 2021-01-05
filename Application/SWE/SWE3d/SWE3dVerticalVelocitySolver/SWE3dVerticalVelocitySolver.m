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
        % This corresponds to the matrix applied to the horizontal partial
        % derivative term
        TempRHSCoeMatrix
        % This corresponds to the matrix applied to the vertical velocity
        % term
        TempVertCoeMatrix
        StiffMatrix
        InnerEdge2d
        BoundaryEdge2d
        InnerEdge3d
        BoundaryEdge3d
        mesh2d
        mesh3d
        cell2d
        cell3d
    end
    
    methods
        function obj = SWE3dVerticalVelocitySolver( mesh2d, mesh3d)
            
            warning('off');
            obj.InnerEdge2d = struct(mesh2d.InnerEdge);
            obj.BoundaryEdge2d = struct(mesh2d.BoundaryEdge);
            obj.InnerEdge3d = struct(mesh3d.InnerEdge);
            obj.BoundaryEdge3d = struct(mesh3d.BoundaryEdge);
            obj.mesh2d = struct(mesh2d);
            obj.mesh3d = struct(mesh3d);
            obj.cell2d = struct(mesh2d.cell);
            obj.cell3d = struct(mesh3d.cell);
            warning('on');
            
            obj.matClearGlobalMemory( );
            
            obj.RHSCoeMatrix = cell(1);
            obj.RHSCoeMatrix{1} = zeros( mesh3d.cell.Np, mesh3d.cell.Np, mesh3d.K );
            obj.VertCoeMatrix = obj.RHSCoeMatrix;
            %% we calculate from bottom to surface
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
                obj.RHSCoeMatrix{1}(:,:, i*Nz )  = CoeMatrix\(diag(mesh3d.J(:, i*Nz )) * mesh3d.cell.M);
            end
            %% we calculate from surface to bottom
            %             UpEidM = mesh3d.cell.Fmask( mesh3d.cell.Fmask( :,end) ~= 0, end );
            %             Nz = mesh3d.Nz;
            %             for i = 1:mesh2d.K
            %                 M3d = zeros(mesh3d.cell.Np);
            %                 M3d( UpEidM, UpEidM ) = diag(mesh2d.J(:,i)) * mesh2d.cell.M;
            %                 %                             CoeMatrix = diag(mesh3d.J(:, (i-1)*Nz+1)) * mesh3d.cell.M * diag(mesh3d.tz(:, (i-1)*Nz+1)) * mesh3d.cell.Dt -  M3d;
            %                 % if the outer value at the surface is taken to be the opposite value of the inner data, then the CoeMatrix for the uppermost cell is set to be the following
            %                 CoeMatrix = diag(mesh3d.J(:, (i-1)*Nz+1)) * mesh3d.cell.M * diag(mesh3d.tz(:, (i-1)*Nz+1)) * mesh3d.cell.Dt -  2*M3d;
            %                 obj.RHSCoeMatrix{1}(:,:, (i-1)*Nz+1 )  = CoeMatrix\(diag(mesh3d.J(:, (i-1)*Nz+1 )) * mesh3d.cell.M);
            %                 %for the uppermost cell, the numerical flux for omega is zero, so the VertCoeMatrix for this cell is set to be zero
            %                 for j = 2 : Nz
            %                     CoeMatrix = diag(mesh3d.J(:, (i-1)*Nz + j)) * mesh3d.cell.M * diag(mesh3d.tz(:, (i-1)*Nz + j)) * mesh3d.cell.Dt - M3d;
            %                     obj.RHSCoeMatrix{1}(:,:,(i-1)*Nz + j)  = CoeMatrix\(diag(mesh3d.J(:, (i-1)*Nz + j)) * mesh3d.cell.M);
            %                     obj.VertCoeMatrix{1}(:,:,(i-1)*Nz + j)  = CoeMatrix\(-1*M3d);
            %                 end
            %             end
            %% we calculate the vertical velocity globally
%                         obj.StiffMatrix = cell(1);
%                         obj.StiffMatrix{1} = zeros(mesh3d.cell.Np*mesh3d.Nz, mesh3d.cell.Np*mesh3d.Nz, mesh2d.K);
%                         BottomEidM   = mesh3d.cell.Fmask(mesh3d.cell.Fmask(:,end-1)~=0,end-1);
%                         UpEidM     = mesh3d.cell.Fmask(mesh3d.cell.Fmask(:,end)~=0,end);
%                         Np = mesh3d.cell.Np;
%                         Nz = mesh3d.Nz;
%                         for i =1:mesh2d(1).K
%                             %> At present, we assume the mesh is uniform in the vertical direction
%                             ElementalMassMatrix3d = diag(mesh3d.J(:,(i-1)*Nz+1)) * mesh3d.cell.M;
%                             ElementalMassMatrix2d = diag(mesh2d.J(:,i)) * mesh2d.cell.M;
%                             Dz3d = diag(mesh3d.tz(:,(i-1)*Nz+1)) * mesh3d.cell.Dt;
%                             %> first cell first, then the other cells left
%                             LocalRows    = (1:Np)';
%                             LocalColumns = 1:Np;
%                             OP11 = ElementalMassMatrix3d * Dz3d;
%                             OP11 = ImposeSurfaceDirichletBoundary(OP11, UpEidM, ElementalMassMatrix2d);
%                             if mesh3d.Nz ~= 1
%                                 OP11 = LocalDownBoundaryIntegral(BottomEidM, ElementalMassMatrix2d, OP11);
%                                 obj.StiffMatrix{1}(LocalRows,LocalColumns,i) = ElementalMassMatrix3d\OP11;
%                                 AdjacentRows = (Np+1:2*Np)';
%                                 OP12 = zeros(Np);
%                                 OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, OP12, ElementalMassMatrix2d);
%                                 obj.StiffMatrix{1}(AdjacentRows,LocalColumns,i) = ElementalMassMatrix3d\OP12;
%                                 for j = 2:mesh3d.Nz-1
%                                     UpAdjacentRows = ((j-2)*Np+1:(j-1)*Np)';
%                                     LocalRows    = ((j-1)*Np+1:j*Np)';
%                                     BottomAdjacentRows = (j*Np+1:(j+1)*Np)';
%                                     LocalColumns = (j-1)*Np+1:j*Np;
%                                     OP11 = ElementalMassMatrix3d * Dz3d;
%                                     OP11 = LocalUpBoundaryIntegral(UpEidM, ElementalMassMatrix2d, OP11);
%                                     OP11 = LocalDownBoundaryIntegral(BottomEidM, ElementalMassMatrix2d, OP11);
%                                     obj.StiffMatrix{1}(LocalRows,LocalColumns,i) = ElementalMassMatrix3d\OP11;
%                                     OP12 = zeros(Np);
%                                     OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, OP12, ElementalMassMatrix2d);
%                                     obj.StiffMatrix{1}(UpAdjacentRows,LocalColumns,i) = ElementalMassMatrix3d\OP12;
%                                     OP12 = zeros(Np);
%                                     OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, OP12, ElementalMassMatrix2d);
%                                     obj.StiffMatrix{1}(BottomAdjacentRows,LocalColumns,i) = ElementalMassMatrix3d\OP12;
%                                 end
%                                 LocalRows    = ((Nz-1)*Np+1:Nz*Np)';
%                                 LocalColumns = ((Nz-1)*Np+1:Nz*Np);
%                                 AdjacentRows = ((Nz-2)*Np+1:(Nz-1)*Np)';
%                                 OP11 = ElementalMassMatrix3d * Dz3d;
%                                 OP11 = ImposeBottomDirichletBoundary(OP11, BottomEidM, ElementalMassMatrix2d);
%                                 obj.StiffMatrix{1}(LocalRows,LocalColumns,i) = ElementalMassMatrix3d\OP11;
%                                 OP12 = zeros(Np);
%                                 OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, OP12, ElementalMassMatrix2d);
%                                 obj.StiffMatrix{1}(AdjacentRows,LocalColumns,i) = ElementalMassMatrix3d\OP12;
%                             else
%                                 OP11 = ImposeBottomDirichletBoundary(OP11, BottomEidM, ElementalMassMatrix2d);
%                                 obj.StiffMatrix{1}(LocalRows,LocalColumns,i) = ElementalMassMatrix3d\OP11;
%                             end
%                         end
        end
        
        function VerticalVelocity = matCalculateVerticalVelocity( obj, physClass, fphys2d, fphys )
            BotEidM = physClass.meshUnion.cell.Fmask( physClass.meshUnion.cell.Fmask( :,end-1) ~= 0, end-1 );
            UpEidM = physClass.meshUnion.cell.Fmask( physClass.meshUnion.cell.Fmask( :,end) ~= 0, end );  
            % The C version from bottom to surface      
% %             tic;
            [ VerticalVelocity ] = mxCalculateVerticalVelocity(obj.mesh2d, obj.mesh3d, obj.InnerEdge2d, obj.BoundaryEdge2d, obj.InnerEdge3d,  obj.BoundaryEdge3d,...
                fphys2d{1}, fphys{1}, physClass.hcrit, obj.cell2d, obj.cell3d, physClass.gra, physClass.fext2d{1}, physClass.fext3d{1}, obj.RHSCoeMatrix{1},...
                obj.VertCoeMatrix{1}, BotEidM, UpEidM, int8(physClass.meshUnion.mesh2d.BoundaryEdge.ftype), int8(physClass.meshUnion.BoundaryEdge.ftype) );
% %             t1 = toc;
% %             tic;
% %             edge = physClass.meshUnion.InnerEdge;
% %             edge2d = physClass.mesh2d.InnerEdge;
% %             [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
% %             [ FluxM3d ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
% %             %             [ FluxM2d ] = edge.VerticalColumnIntegralField(FluxM3d(:,:,1) );
% %             [ FluxP3d ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
% %             %             [ FluxP2d ] = edge.VerticalColumnIntegralField(FluxP3d(:,:,1) );
% %             [ FluxS3d ] = physClass.matEvaluateSurfNumFlux( physClass.meshUnion, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), fp(:,:,[4, 1, 2]), edge );
% %             %             [ FluxS2d ] = edge.VerticalColumnIntegralField(FluxS3d(:,:,1) );
% %             [ InnerSurface_frhs3d ] = edge.matEvaluateStrongFromEdgeRHS( FluxM3d(:,:,1), FluxP3d(:,:,1), FluxS3d(:,:,1) ) ;            
% %             [ fm2d, fp2d ] = edge2d.matEvaluateSurfValue( fphys2d );
% %             [ fluxM ] = physClass.PCESolver2d.matEvaluateSurfFlux( edge2d, edge2d.nx, edge2d.ny, fm2d );
% %             [ fluxP ] = physClass.PCESolver2d.matEvaluateSurfFlux( edge2d, edge2d.nx, edge2d.ny, fp2d );
% %             [ fluxS ] = physClass.matEvaluateSurfNumFlux( physClass.mesh2d, edge2d.nx, edge2d.ny, fm2d, fp2d, edge2d);
% %             InnerSurface_frhs2d = edge2d.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS(:,:,1) );
% %             
% %             %             [ InnerSurface_frhs2d ] = edge2d.matEvaluateStrongFormEdgeRHS( FluxM2d(:,:,1), FluxP2d(:,:,1), FluxS2d(:,:,1) ) ;
% %             
% %             edge = physClass.meshUnion.BoundaryEdge;
% %             edge2d = physClass.mesh2d.BoundaryEdge;
% %             [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
% %             [ fm, fp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, fm, fp, physClass.fext3d{1} );
% %             [ FluxM3d ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
% %             %             [ FluxM2d ] = edge.VerticalColumnIntegralField( FluxM3d(:,:,1) );
% %             [ FluxS3d ] = physClass.matEvaluateSurfNumFlux( physClass.meshUnion, edge.nx, edge.ny, fm(:,:,[4, 1, 2]), fp(:,:,[4, 1, 2]), edge );
% %             %             [ FluxS2d ] = edge.VerticalColumnIntegralField( FluxS3d(:,:,1) );
% %             
% %             [ BoundarySurface_frhs3d ] = edge.matEvaluateStrongFormEdgeRHS( FluxM3d(:,:,1), FluxS3d(:,:,1) ) ;
% % %             [ BoundarySurface_frhs2d ] = edge2d.matEvaluateStrongFormEdgeRHS( FluxM2d(:,:,1), FluxS2d(:,:,1) );
% %             [ fm2d, fp2d ] = edge2d.matEvaluateSurfValue( fphys2d );
% %             [ fm2d, fp2d ] = physClass.PCESolver2d.matImposeBoundaryCondition( edge2d, edge2d.nx, edge2d.ny, fm2d, fp2d, physClass.fext2d{1});            
% %             [ fluxM ] = physClass.PCESolver2d.matEvaluateSurfFlux( edge2d, edge2d.nx, edge2d.ny, fm2d );
% %             [ fluxS ] = physClass.matEvaluateSurfNumFlux( physClass.mesh2d, edge2d.nx, edge2d.ny, fm2d, fp2d, edge2d);
% %             BoundarySurface_frhs2d = edge2d.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
% %             
% %             field3d =...
% %                 ( - physClass.meshUnion.rx .* (physClass.meshUnion.cell.Dr * fphys{1}(:,:,1)) - physClass.meshUnion.sx .* (physClass.meshUnion.cell.Ds * fphys{1}(:,:,1)) ) + ...
% %                 InnerSurface_frhs3d + BoundarySurface_frhs3d  - ( physClass.meshUnion.ry .* (physClass.meshUnion.cell.Dr * fphys{1}(:,:,2)) - physClass.meshUnion.sy .* (physClass.meshUnion.cell.Ds * fphys{1}(:,:,2)) ) ;
% %             % Term2d = ;
% %             field2d = physClass.meshUnion.Extend2dField( - physClass.mesh2d.rx .* (physClass.mesh2d.cell.Dr * fphys2d{1}(:,:,2)) - physClass.mesh2d.sx .* (physClass.mesh2d.cell.Ds * fphys2d{1}(:,:,2) ) + ...
% %                 InnerSurface_frhs2d + BoundarySurface_frhs2d - physClass.mesh2d.ry .* (physClass.mesh2d.cell.Dr * fphys2d{1}(:,:,3) ) - physClass.mesh2d.sy .* (physClass.mesh2d.cell.Ds * fphys2d{1}(:,:,3)) );
% %             
% %             %             x = physClass.meshUnion.x;
% % %             y = physClass.meshUnion.y;
% % %             z = physClass.meshUnion.z;
% % %             t = time;
% % %             field2d = field2d  + eval( physClass.Source2d );  
% %             
% %             Nz = physClass.meshUnion.Nz;
% %             VerticalVelocity = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
% %             %% we calculate from bottom to surface
% %             VerticalVelocity(:, Nz:Nz:end) =  permute( sum( bsxfun(@times, obj.RHSCoeMatrix{1}(:,:,Nz:Nz:end), ...
% %                 permute( permute( field3d(:, Nz:Nz:end) - field2d(:, Nz:Nz:end), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
% %             
% %             for Layer = 2 : Nz
% %                 BotVertVelocity = zeros( physClass.meshUnion.cell.Np, physClass.mesh2d.K );
% %                 BotVertVelocity( BotEidM, : ) = VerticalVelocity( UpEidM, Nz - Layer + 2 : Nz : end);
% %                 VerticalVelocity(:, Nz - Layer + 1:Nz:end) = permute( sum( bsxfun(@times, obj.RHSCoeMatrix{1}(:,:,Nz - Layer + 1:Nz:end), ...
% %                     permute( permute( field3d(:, Nz - Layer + 1:Nz:end) - field2d(:, Nz - Layer + 1:Nz:end), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]) + permute( sum( bsxfun(@times, obj.VertCoeMatrix{1}(:,:,Nz - Layer + 1:Nz:end), ...
% %                     permute( permute( BotVertVelocity, [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
% %             end
% %             
% %             t2 = toc;
% %             fprintf("The speed up ratio is:%f\n", t2/t1);
            
            %% we calculate from surface to bottom
            %             VerticalVelocity = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
            %             VerticalVelocity(:, 1:Nz:end) =  permute( sum( bsxfun(@times, obj.RHSCoeMatrix{1}(:,:,1:Nz:end), ...
            %                 permute( permute( field3d(:, 1:Nz:end) - field2d(:, 1:Nz:end), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
            %
            %             for Layer = 2 : Nz
            %                 UpVertVelocity = zeros( physClass.meshUnion.cell.Np, physClass.mesh2d.K );
            %                 UpVertVelocity( UpEidM, : ) = VerticalVelocity( BotEidM, Layer - 1 : Nz : end);
            %                 VerticalVelocity(:, Layer:Nz:end) = permute( sum( bsxfun(@times, obj.RHSCoeMatrix{1}(:,:,Layer:Nz:end), ...
            %                     permute( permute( field3d(:, Layer:Nz:end) - field2d(:, Layer:Nz:end), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]) + permute( sum( bsxfun(@times, obj.VertCoeMatrix{1}(:,:,Layer:Nz:end), ...
            %                     permute( permute( UpVertVelocity, [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
            %             end
            %% The vertical velocity is calculated implicitly with central flux adopted
            %             VerticalVelocity = zeros( physClass.meshUnion.cell.Np, physClass.meshUnion.K );
            %             for i = 1:physClass.mesh2d.K
            %                 tempdata = field3d(:, (i-1)*Nz+1:i*Nz) - field2d(:, (i-1)*Nz+1:i*Nz);
            %                 VerticalVelocity((i-1)*Nz*physClass.meshUnion.cell.Np+1:i*Nz*physClass.meshUnion.cell.Np) = obj.StiffMatrix{1}(:,:,i)\tempdata(:);
            %             end
        end
        
        function matClearGlobalMemory(obj)
            clear mxCalculateVerticalVelocity;
        end
    end
    
end

function OP11 = ImposeSurfaceDirichletBoundary(OP11, UpEidM, ElementalMassMatrix2d)
OP11(UpEidM, UpEidM) = OP11(UpEidM, UpEidM) - ElementalMassMatrix2d;
end

function OP11 = ImposeBottomDirichletBoundary(OP11, BottomEidM, ElementalMassMatrix2d)
OP11(BottomEidM, BottomEidM) = OP11(BottomEidM, BottomEidM) + ElementalMassMatrix2d;
end

function OP11 = LocalDownBoundaryIntegral(BottomEidM, ElementalMassMatrix2d, OP11)
OP11(BottomEidM, BottomEidM) = OP11(BottomEidM, BottomEidM) + 1/2*ElementalMassMatrix2d;
end

function OP11 = LocalUpBoundaryIntegral(UpEidM, ElementalMassMatrix2d, OP11)
OP11(UpEidM,UpEidM) = OP11(UpEidM,UpEidM) - 1/2*ElementalMassMatrix2d;
end

function OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, OP12, ElementalMassMatrix2d)
OP12(BottomEidM, UpEidM) = OP12(BottomEidM, UpEidM) + 1/2*ElementalMassMatrix2d;
end

function OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, OP12, ElementalMassMatrix2d)
OP12(UpEidM, BottomEidM) = OP12(UpEidM, BottomEidM) - 1/2*ElementalMassMatrix2d;
end


