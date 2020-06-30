classdef NdgHorizDiffSolver < AbstractDiffSolver
    %NDGHORIZDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties( Access = protected )
        px
        py
        InnerEdgeTau
        BoundaryEdgeTau
        M
        invM
    end
    
    methods
        function obj = NdgHorizDiffSolver( physClass )
            obj = obj@AbstractDiffSolver( physClass );
            obj.assembleMassMatrix( physClass.meshUnion(1) );
            if physClass.option.isKey('AdvDiffHorizontalDiffusionType')
               if physClass.option.isKey('AdvDiffConstantHorizontalDiffusionValue')
                   value = physClass.getOption('AdvDiffConstantHorizontalDiffusionValue');
                   obj.nv = value * ones(size(physClass.meshUnion(1).x));
                   fprintf('Value of the constant horizontal diffusion coefficient is set to be: %f\n',value);
                   obj.matUpdatePenaltyParameter(  physClass, obj.nv );
               end
            end            
        end
        
        function matEvaluateDiffRHS(obj, physClass, fphys)
            Kappa = obj.nv;
            for i = 1:physClass.Nvar
                obj.matCalculateAuxialaryVariable( physClass, fphys(:,:,physClass.varFieldIndex(i)), Kappa, i, ...
                    physClass.InnerEdgefm{1}(:,:,physClass.varFieldIndex(i)),...
                    physClass.InnerEdgefp{1}(:,:,physClass.varFieldIndex(i)), ...
                    physClass.BoundaryEdgefm{1}(:,:,physClass.varFieldIndex(i)),...
                    physClass.BoundaryEdgefp{1}(:,:,physClass.varFieldIndex(i)));
            end
            %> this part is used to calculate $\frac{\partial}{\partial x}(\nv_h \frac{\partial c}{\partial x})
            %> + \frac{\partial}{\partial y}(\nv (\frac{\partial c}{\partial y}))$
            for i = 1:physClass.Nvar
            physClass.frhs{1}(:,:,i) = physClass.frhs{1}(:,:,i) + obj.matCalculatePartDerivTermX( physClass, obj.px(:,:,i),...
                Kappa, fphys(:,:,physClass.varFieldIndex(i)), 1, physClass.InnerEdgefm{1}(:,:,physClass.varFieldIndex(i)), ...
                physClass.InnerEdgefp{1}(:,:,physClass.varFieldIndex(i)), physClass.BoundaryEdgefm{1}(:,:,physClass.varFieldIndex(i)),...
                physClass.BoundaryEdgefp{1}(:,:,physClass.varFieldIndex(i))) + ...
                obj.matCalculatePartDerivTermY( physClass, obj.py(:,:,i),...
                Kappa, fphys(:,:,physClass.varFieldIndex(i)), 1, physClass.InnerEdgefm{1}(:,:,physClass.varFieldIndex(i)), ...
                physClass.InnerEdgefp{1}(:,:,physClass.varFieldIndex(i)), physClass.BoundaryEdgefm{1}(:,:,physClass.varFieldIndex(i)),...
                physClass.BoundaryEdgefp{1}(:,:,physClass.varFieldIndex(i)));
           
            end
       
        end        
    end
    
    methods( Access = protected )
        
        function matUpdateViscosity(obj)
            %doing nothing
        end
        
        function [ invM ] = assembleMassMatrix( obj, mesh )
            cell = mesh.cell;
            Np = cell.Np;
            K = mesh.K;
            invM = zeros( Np, Np, K );
            for k = 1:K
                Jq = cell.project_node2quad( mesh.J(:, k) );
                obj.M(:,:,k) = ( cell.Vq' * diag( Jq.*cell.wq ) ) * cell.Vq;
                obj.invM(:, :, k) = inv( obj.M(:,:,k) );
            end
        end
        
        function matCalculateAuxialaryVariable(obj, physClass, fphys, Kappa, VarIndex, InnerEdgefm, InnerEdgefp, BoundaryEdgefm, BoundaryEdgefp)
            %> @brief Calculating the auxialary variable, i.e. the first order partial derivative of variable U
            %> @detail this function is used to calculate the auxialary variable corresponding to
            %> the first order partial derivative  of variable U, i.e. it calculates
            %> $q_x=\Kappa\frac{\partial U}{\partial x}$ and $q_y=\Kappa\frac{\partial U}{\partial y}$
            %> @param[in] physClass The physical solver establised
            %> @param[in] fphys The physical field to be calculated, i.e. U in the presentation
            %> @param[in] Kappa The diffusion coefficient
            %> @param[in] VarIndex The index of the calculated variable, it
            %> determines the order of $q_x$ and $q_y$ to be placed in obj.px
            %> and obj.py.
            %> @param[in] InnerEdgefm The local value of fphys at InnerEdge
            %> @param[in] InnerEdgefp The adjacent value of fphys at InnerEdge
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge,
            %> we point out that both the local value and the adjacent value at
            %> the boundary edge and inner edge are input from the physical solver
            %> establised, this is because these variables are frequently
            %> used in the current three-dimensional solver, and we don't
            %> want to fetch it every time when needed and also don't want to impose
            %> the boundary condition, so we choose to fetch once and store it
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
%                         for k = 1:physClass.meshUnion(1).K
%                             % The diffusion coeffcient has been grouped into the test
%                             % function, which means the massmatrix need to multiply the
%                             % diffusion coefficient
%                             obj.px(:,k,VarIndex) = obj.invM(:,:,k) * diag(Kappa(:,k)) *...
%                                 obj.M(:,:,k) * ( physClass.meshUnion(1).rx(:,k) .* (physClass.meshUnion(1).cell.Dr *  fphys(:,k) ) + ...
%                             physClass.meshUnion(1).sx(:,k) .* ( physClass.meshUnion(1).cell.Ds *  fphys(:,k) ));
%                             obj.py(:,k,VarIndex) = obj.invM(:,:,k) * diag(Kappa(:,k)) *...
%                                 obj.M(:,:,k) * ( physClass.meshUnion(1).ry(:,k) .* (physClass.meshUnion(1).cell.Dr *  fphys(:,k) ) + ...
%                             physClass.meshUnion(1).sy(:,k) .* ( physClass.meshUnion(1).cell.Ds *  fphys(:,k) ));
%                         end

            obj.px(:,:,VarIndex) =  permute( sum( bsxfun(@times, obj.invM, ...
                permute( permute( Kappa .* permute( sum( bsxfun(@times, obj.M, ...
                permute( permute( ( physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr *  fphys ) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds *  fphys )), [1,3,2] ), ...
                [2,1,3] ) ), 2 ), [1,3,2]), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
            
            obj.py(:,:,VarIndex) = permute( sum( bsxfun(@times, obj.invM, ...
                permute( permute( Kappa .* permute( sum( bsxfun(@times, obj.M, ...
                permute( permute( ( physClass.meshUnion(1).ry .* (physClass.meshUnion(1).cell.Dr *  fphys ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds *  fphys )), [1,3,2] ), ...
                [2,1,3] ) ), 2 ), [1,3,2]), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
            edge = physClass.meshUnion(1).InnerEdge;
            %             [fM, fP] = obj.matEvaluateSurfValue(edge, fphys );
            [KappaM, KappaP] = obj.matEvaluateSurfValue(edge, Kappa );
%                             [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( physClass.InnerEdgeFluxM3d{m}(:,:,[2,3]), physClass.InnerEdgeFluxP3d{m}(:,:,[2,3]), physClass.InnerEdgeFluxS3d{m}(:,:,[2,3]) );

            obj.px(:,:,VarIndex) = obj.px(:,:,VarIndex) - edge.matEvaluateStrongFromEdgeRHS( ...
                KappaM .* InnerEdgefm .* edge.nx, KappaP .* InnerEdgefp .* edge.nx , 0.5 * (KappaM .* InnerEdgefm + KappaP .* InnerEdgefp) .* edge.nx);
            obj.py(:,:,VarIndex) = obj.py(:,:,VarIndex) - edge.matEvaluateStrongFromEdgeRHS( ...
                KappaM .* InnerEdgefm .* edge.ny, KappaP .* InnerEdgefp .* edge.ny , 0.5 * (KappaM .* InnerEdgefm + KappaP .* InnerEdgefp) .* edge.ny);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            %             [fM, fP] = obj.matEvaluateSurfValue(edge, fphys );
            [KappaM, ~] = obj.matEvaluateSurfValue(edge, Kappa );
            obj.px(:,:,VarIndex) = obj.px(:,:,VarIndex) - edge.matEvaluateStrongFormEdgeRHS( ...
                KappaM .* BoundaryEdgefm .* edge.nx , 0.5 * KappaM .* (BoundaryEdgefm + BoundaryEdgefp) .* edge.nx );
            obj.py(:,:,VarIndex) = obj.py(:,:,VarIndex) - edge.matEvaluateStrongFormEdgeRHS( ...
                KappaM .* BoundaryEdgefm .* edge.ny , 0.5 * KappaM .* (BoundaryEdgefm + BoundaryEdgefp) .* edge.ny);
        end
        
        function frhs = matCalculatePartDerivTermX(obj, physClass, pfield, Kappa, fphys, Prantl, InnerEdgefm, InnerEdgefp, BoundaryEdgefm, BoundaryEdgefp)
            %> @brief Calculating the second order derivative in x direction
            %> @detail this function is used to calculate the second order
            %> derivative in x direction, i.e. $\frac{\partial \sigma}{\partial x}$,
            %> in the above formulation, $\sigma$ can be $\Kappa\frac{\partial U}{\partial x}$
            %> or $\Kappa\frac{\partial U}{\partial y}$
            %> @param[in] physClass The physical solver establised
            %> @param[in] pfield Value of the first order derivative of U,
            %> this value has been contained in obj.px or obj.py, we just
            %> need to fetch it as we need
            %> @param[in] Kappa The diffusion coefficient, this is needed to
            %> calculate the numerical flux
            %> @param[in] fphys The physical field to be calculated, i.e. U
            %> in the presentation
            %> @param[in] Prantl The Prantl number, for conventional convection-diffusion
            %> problem, this value equals to 1, however, for salt and tempreture contained
            %> in SWE3d, this value equals to 1.1, this value is needed here, since it
            %> influence Kappa and the penalty parameter.
            %> @param[in] InnerEdgefp The local value of fphys at InnerEdge, used to
            %> calculate numerical flux
            %> @param[in] InnerEdgefp The adjacent value of fphys at InnerEdge, used
            %> to calculate numerical flux
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge,
            %> we point out that both the local value and the adjacent value at
            %> the boundary edge and inner edge are input from the physical solver
            %> establised, this is because these variables are frequently
            %> used in the current three-dimensional solver, and we don't
            %> want to fetch it every time when needed and also don't want to impose
            %> the boundary condition, so we choose to fetch once and store it
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
            
            frhs = physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr * pfield ) + ...
                physClass.meshUnion(1).sx .* (physClass.meshUnion(1).cell.Ds * pfield );
            
%             LocalVariable = zeros( physClass.meshUnion(1).cell.Np, physClass.meshUnion(1).K );
%             for k = 1:physClass.meshUnion(1).K
%                 LocalVariable(:,k) = obj.invM(:,:,k) * diag(Kappa(:,k)) *...
%                     obj.M(:,:,k) * ( physClass.meshUnion(1).rx(:,k) .* (physClass.meshUnion(1).cell.Dr *  fphys(:,k) ) + ...
%                     physClass.meshUnion(1).sx(:,k) .* ( physClass.meshUnion(1).cell.Ds *  fphys(:,k) ));
%             end
            LocalVariable = permute( sum( bsxfun(@times, obj.invM, ...
                permute( permute( Kappa .* permute( sum( bsxfun(@times, obj.M, ...
                permute( permute( ( physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr *  fphys ) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds *  fphys )), [1,3,2] ), ...
                [2,1,3] ) ), 2 ), [1,3,2]), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);
            
            %             LocalVariable = Kappa .* ( physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr * fphys ) + ...
            %                 physClass.meshUnion(1).sx .* (physClass.meshUnion(1).cell.Ds * fphys ));
            
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ fluxP ] = pfieldP .* edge.nx;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.nx .* (LocalVariableM + LocalVariableP)./2 - edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.nx - InnerEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            %             [ fluxP ] = pfieldP .* edge.nx;
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.nx .* LocalVariableM - edge.nx .* obj.BoundaryEdgeTau./Prantl .* ( BoundaryEdgefm .* edge.nx - BoundaryEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM,  fluxS);
        end
        
        function frhs = matCalculatePartDerivTermY(obj, physClass, pfield, Kappa, fphys, Prantl, InnerEdgefm, InnerEdgefp, BoundaryEdgefm, BoundaryEdgefp)
            %> @brief Calculating the second order derivative in y direction
            %> @detail this function is used to calculate the second order
            %> derivative in y direction, i.e. $\frac{\partial \sigma}{\partial y}$,
            %> in the above formulation, $\sigma$ can be $\Kappa\frac{\partial U}{\partial y}$
            %> or $\Kappa\frac{\partial U}{\partial x}$
            %> @param[in] physClass The physical solver establised
            %> @param[in] pfield Value of the first order derivative of U,
            %> this value has been contained in obj.px or obj.py, we just
            %> need to fetch it as we need
            %> @param[in] Kappa The diffusion coefficient, this is needed to
            %> calculate the numerical flux
            %> @param[in] fphys The physical field to be calculated, i.e. U
            %> in the presentation
            %> @param[in] Prantl The Prantl number, for conventional convection-diffusion
            %> problem, this value equals to 1, however, for salt and tempreture contained
            %> in SWE3d, this value equals to 1.1, this value is needed here, since it
            %> influence Kappa and the penalty parameter.
            %> @param[in] InnerEdgefp The local value of fphys at InnerEdge, used to
            %> calculate numerical flux
            %> @param[in] InnerEdgefp The adjacent value of fphys at InnerEdge, used
            %> to calculate numerical flux
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge,
            %> we point out that both the local value and the adjacent value at
            %> the boundary edge and inner edge are input from the physical solver
            %> establised, this is because these variables are frequently
            %> used in the current three-dimensional solver, and we don't
            %> want to fetch it every time when needed and also don't want to impose
            %> the boundary condition, so we choose to fetch once and store it
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
            frhs = physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * pfield ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * pfield );
            
%             LocalVariable = zeros( physClass.meshUnion(1).cell.Np, physClass.meshUnion(1).K );
%             for k = 1:physClass.meshUnion(1).K
%                 LocalVariable(:,k) = obj.invM(:,:,k) * diag(Kappa(:,k)) *...
%                     obj.M(:,:,k) * ( physClass.meshUnion(1).ry(:,k) .* (physClass.meshUnion(1).cell.Dr *  fphys(:,k) ) + ...
%                     physClass.meshUnion(1).sy(:,k) .* ( physClass.meshUnion(1).cell.Ds *  fphys(:,k) ));
%             end
            
            LocalVariable = permute( sum( bsxfun(@times, obj.invM, ...
                permute( permute( Kappa .* permute( sum( bsxfun(@times, obj.M, ...
                permute( permute( ( physClass.meshUnion(1).ry .* (physClass.meshUnion(1).cell.Dr *  fphys ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds *  fphys )), [1,3,2] ), ...
                [2,1,3] ) ), 2 ), [1,3,2]), [1,3,2] ), [2,1,3] ) ), 2 ), [1,3,2]);            
            
            %             LocalVariable = Kappa .* ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * fphys ) + ...
            %                 physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * fphys ) );
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ fluxP ] = pfieldP .* edge.ny;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.ny .* (LocalVariableM + LocalVariableP)./2 - ...
                edge.ny .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.ny - InnerEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            %             [ fluxP ] = pfieldP .* edge.ny;
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.ny .* LocalVariableM  - edge.ny .* obj.BoundaryEdgeTau./Prantl .* ( BoundaryEdgefm .* edge.ny - BoundaryEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS);
        end
        
        function matUpdatePenaltyParameter( obj, physClass, DiffusionCoefficient )
            %> this penalty parameter is calculated as $\tau=\frac{(D_p+1)(D_p+d)}{d}\frac{n_0}{2}\frac{A}{V}\miu$
            [ HnvM, HnvP ] = obj.matEvaluateSurfValue(physClass.meshUnion(1).InnerEdge, DiffusionCoefficient );
            
            InnerEdgeA_fm = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(1,:)))', 1, physClass.meshUnion(1).Nz );
            InnerEdgeA_fp = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(2,:)))', 1, physClass.meshUnion(1).Nz );
            InnerEdgeTau_fm = bsxfun(@times,  ( InnerEdgeA_fm(:) )',...
                ( physClass.meshUnion(1).cell.N + 1 )*(physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )/double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* HnvM);
            InnerEdgeTau_fp = bsxfun(@times,  ( InnerEdgeA_fp(:) )',...
                ( physClass.meshUnion(1).cell.N + 1 )*(physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )/double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* HnvP);
            obj.InnerEdgeTau = max( InnerEdgeTau_fm, InnerEdgeTau_fp );
            
            BoundaryEdgeA_fm = repmat ( (physClass.mesh2d(1).BoundaryEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).BoundaryEdge.FToE(1,:)))', 1, physClass.meshUnion(1).Nz );
            [ Hnv, ~ ] = obj.matEvaluateSurfValue(physClass.meshUnion(1).BoundaryEdge, DiffusionCoefficient);
            obj.BoundaryEdgeTau = bsxfun(@times, ( BoundaryEdgeA_fm(:) )', ...
                ( physClass.meshUnion(1).cell.N + 1 )*( physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )./double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* Hnv);
        end        
    end
end

