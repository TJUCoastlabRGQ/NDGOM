classdef NdgQuadFreeStrongLDGVisSolver2d < NdgAbstractVisSolver
    %NDGQUADFREESTRONGLDGVISSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        tau = 1.0
        Beta = [1.0 1.0]
    end
    
    methods
        function obj = NdgQuadFreeStrongLDGVisSolver2d( phys, varId, rhsId )
            obj = obj@NdgAbstractVisSolver( phys, varId, rhsId );
        end
        
        function matEvaluateRHS( obj, fphys )
            matEvaluateAuxiVar( obj, fphys );
            matEvaluateOriVarRHS( obj, fphys );
        end
    end
    
    methods ( Access = protected )
        function matEvaluateAuxiVar( obj, fphys )
            matEvaluateAuxiVarVolumeKernel( obj, fphys );
            matEvaluateAuxiVarSurfaceKernel( obj, fphys );
            % todo: consider boundary condition
        end
        
        function [varx, vary] = matEvaluateScalarJump( obj, Edge, fphys )
             for fld = obj.Nfield
                 id = obj.varId(fld);
                 [ fm, fp ] = Edge.matEvaluateSurfValue( fphys );
                 varx = Edge.nx .* fm(:,:,id) - Edge.nx .* fp(:,:,id);
                 vary = Edge.ny .* fm(:,:,id) - Edge.ny .* fp(:,:,id);
             end
        end
        
        function JumpTerm = matEvaluateVectorJump( obj, Edge )
            [ pxm, pxp ] = Edge.matEvaluateSurfValue( obj.px );
            [ pym, pyp ] = Edge.matEvaluateSurfValue( obj.py );
            JumpTerm = Edge.nx .* pxm - Edge.nx .* pxp + ...
                Edge.ny .* pym - Edge.ny .* pyp;
        end
        
        function matEvaluateAuxiVarVolumeKernel( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.phys.meshUnion(m);
                for fld = obj.Nfield
                    id = obj.varId(fld);
                    dfdr = mesh.cell.Dr * fphys{m}(:, :, id);
                    dfds = mesh.cell.Ds * fphys{m}(:, :, id);
                    
                    obj.px{m}(:, :, fld) = - obj.mx{m} .* ( ...
                        mesh.rx .* dfdr + mesh.sx .* dfds );
                    
                    obj.py{m}(:, :, fld) = - obj.my{m} .* ( ...
                        mesh.ry .* dfdr + mesh.sy .* dfds );
                end
            end
        end
        
        function matEvaluateAuxiVarSurfaceKernel( obj, fphys )
            for m = 1:obj.Nmesh
                edge = obj.phys.meshUnion(m).InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                
                fluxM = fm(:, :, obj.varId) .* edge.nx;
                fluxP = fp(:, :, obj.varId) .* edge.nx;
                % fluxS = ( fluxM + fluxP ) * 0.5;
%                 [Jumpx, Jumpy] = obj.matEvaluateScalarJump( edge, fphys );
%                 fluxS = (fluxM + fluxP) * 0.5  - obj.Beta(1) * edge.nx .* Jumpx - obj.Beta(2) * edge.nx .* Jumpy;
                 fluxS = fm(:, :, obj.varId) .* edge.nx;                

                obj.px{m} = obj.px{m} + obj.mx{m} .* ...
                    edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );
                
                fluxM = fm(:, :, obj.varId) .* edge.ny;
                fluxP = fp(:, :, obj.varId) .* edge.ny;
                % fluxS = ( fluxM + fluxP ) * 0.5;
%                 fluxS = (fluxM + fluxP) * 0.5  - obj.Beta(1) * edge.ny .* Jumpx - obj.Beta(2) * edge.ny .* Jumpy;
                fluxS = fm(:, :, obj.varId) .* edge.ny;        
                
                obj.py{m} = obj.py{m} + obj.my{m} .* ...
                    edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );
            end
        end
        
        function matEvaluateOriVarRHS( obj, fphys )
            matEvaluateOriVarSurfaceKernel( obj, fphys );
            matEvaluateOriVarVolumeKernel( obj );
        end
        
        function matEvaluateOriVarVolumeKernel( obj )
            for m = 1:obj.Nmesh
                mesh = obj.phys.meshUnion(m);
                for fld = 1:obj.Nfield
                    id = obj.rhsId(fld);
                    
                    obj.phys.frhs{m}(:, :, id) = obj.phys.frhs{m}(:, :, id) - ( ...
                        mesh.rx .* ( mesh.cell.Dr * obj.px{m}(:, :, fld) ) + ...
                        mesh.sx .* ( mesh.cell.Ds * obj.px{m}(:, :, fld) ) + ...
                        mesh.ry .* ( mesh.cell.Dr * obj.py{m}(:, :, fld) ) + ...
                        mesh.sy .* ( mesh.cell.Ds * obj.py{m}(:, :, fld) ) );
                end
            end
        end% func
        
        function matEvaluateOriVarSurfaceKernel( obj, fphys )
            for m = 1:obj.Nmesh
                edge = obj.phys.meshUnion(m).InnerEdge;
                [ pxM, pxP ] = edge.matEvaluateSurfValue( obj.px );
                [ pyM, pyP ] = edge.matEvaluateSurfValue( obj.py );
                % [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                % tau = 1;
                for fld = 1:obj.Nfield
                    fluxM = pxM(:, :, fld) .* edge.nx + pyM(:, :, fld) .* edge.ny;
                    fluxP = pxP(:, :, fld) .* edge.nx + pyP(:, :, fld) .* edge.ny;
                end
%                 [Jumpx, Jumpy] = obj.matEvaluateScalarJump( edge, fphys );
%                 VectorJump = obj.matEvaluateVectorJump( edge );
%                 fluxS = ( fluxM + fluxP ) * 0.5 + obj.Beta(1) * edge.nx .* VectorJump + ...
%                     obj.Beta(2) * edge.ny .* VectorJump - obj.tau * edge.nx .* Jumpx - obj.tau * edge.ny .* Jumpy;
                % fluxS = ( fluxM + fluxP ) * 0.5;
                fluxS = pxP(:, :, fld) .* edge.nx + pyP(:, :, fld) .* edge.ny;
                obj.phys.frhs{m}(:, :, obj.rhsId) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId) + ...
                    edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );
            end
        end
    end
    
end

