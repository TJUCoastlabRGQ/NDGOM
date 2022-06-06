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
            %             obj.assembleMassMatrix( physClass.meshUnion(1) );
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
            mesh = physClass.meshUnion(1);
            edge = mesh.InnerEdge;
            [ IEfm, IEfp ] = edge.matEvaluateSurfValue( fphys );
            
            edge = mesh.BoundaryEdge;
            [ BEfm, BEfp ] = edge.matEvaluateSurfValue( fphys );
            
%             [ ~, BEfp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, BEfm, BEfp, physClass.fext{1} );
            [ ~, BEfp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, BEfm, BEfp, physClass.fext{1} );
            Kappa = obj.nv;
            for i = 1:physClass.Nvar
                obj.matCalculateAuxialaryVariable( physClass, fphys{1}(:,:,physClass.varFieldIndex(i)), Kappa, i, ...
                    IEfm(:,:,physClass.varFieldIndex(i)),...
                    IEfp(:,:,physClass.varFieldIndex(i)), ...
                    BEfm(:,:,physClass.varFieldIndex(i)),...
                    BEfp(:,:,physClass.varFieldIndex(i)));
            end
            %> this part is used to calculate $\frac{\partial}{\partial x}(\nu \frac{\partial c}{\partial x})
            %> + \frac{\partial}{\partial y}(\nu (\frac{\partial c}{\partial y}))$
            for i = 1:physClass.Nvar
                physClass.frhs{1}(:,:,i) = physClass.frhs{1}(:,:,i) + obj.matCalculatePartDerivTermX( physClass, obj.px(:,:,i),...
                    Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
                    IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
                    BEfp(:,:,physClass.varFieldIndex(i))) + ...
                    obj.matCalculatePartDerivTermY( physClass, obj.py(:,:,i),...
                    Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
                    IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
                    BEfp(:,:,physClass.varFieldIndex(i)));
                
            end
        end
        
    end
    
    methods( Access = protected )
        
        function matUpdateViscosity(obj)
            %doing nothing
        end
        
%         function [fm, fp] = matImposeBoundaryCondition(obj, physClass, fm, fp)
%             ind = ( edge.ftype == enumBoundaryCondition.Newmann );
%             fp(:, ind, 1) = physClass.GradExt(:, ind, 1);
%             ind = ( edge.ftype == enumBoundaryCondition.Dirichlet );
%             fp(:, ind, 1) = fm(:, ind, 1);
%         end
        
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
            %> $q_x=\kappa\frac{\partial U}{\partial x}$ and $q_y=\kappa\frac{\partial U}{\partial y}$
            %> @param[in] physClass The physical solver establised
            %> @param[in] fphys The physical field to be calculated, i.e. U in the presentation
            %> @param[in] Kappa The diffusion coefficient
            %> @param[in] VarIndex The index of the calculated variable, it
            %> determines the order of $q_x$ and $q_y$ to be placed in obj.px
            %> and obj.py.
            %> @param[in] InnerEdgefm The local value of fphys at InnerEdge
            %> @param[in] InnerEdgefp The adjacent value of fphys at InnerEdge
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
                 
            obj.px(:,:,VarIndex) =  physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr *  fphys ) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds *  fphys );
            
            obj.py(:,:,VarIndex) = physClass.meshUnion(1).ry .* (physClass.meshUnion(1).cell.Dr *  fphys ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds *  fphys );
            
            edge = physClass.meshUnion(1).InnerEdge;
            % The numerical at inner edge is given as $\hat p = (p^-+p^+)/2$
            obj.px(:,:,VarIndex) = obj.px(:,:,VarIndex) - edge.matEvaluateStrongFormEdgeRHS( ...
                InnerEdgefm .* edge.nx, InnerEdgefp .* edge.nx , 0.5 * (InnerEdgefm + InnerEdgefp) .* edge.nx);
            obj.py(:,:,VarIndex) = obj.py(:,:,VarIndex) - edge.matEvaluateStrongFormEdgeRHS( ...
                InnerEdgefm .* edge.ny, InnerEdgefp .* edge.ny , 0.5 * (InnerEdgefm + InnerEdgefp) .* edge.ny);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            % The numerical at Dirichlet boundary is given as $\hat p = p_D$, $\hat p=p$ on Newmann boundary
            % p_D is given when impose boundary condition.
            obj.px(:,:,VarIndex) = obj.px(:,:,VarIndex) - edge.matEvaluateStrongFormEdgeRHS( ...
                BoundaryEdgefm .* edge.nx , BoundaryEdgefp .* edge.nx );
            obj.py(:,:,VarIndex) = obj.py(:,:,VarIndex) - edge.matEvaluateStrongFormEdgeRHS( ...
                BoundaryEdgefm .* edge.ny , BoundaryEdgefp .* edge.ny);
            
            obj.px(:,:,VarIndex) = Kappa .* obj.px(:,:,VarIndex);
            obj.py(:,:,VarIndex) = Kappa .* obj.py(:,:,VarIndex);
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
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
            
            frhs = physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr * pfield ) + ...
                physClass.meshUnion(1).sx .* (physClass.meshUnion(1).cell.Ds * pfield );
            
            LocalVariable = Kappa .* ( physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr * fphys ) + ...
                physClass.meshUnion(1).sx .* (physClass.meshUnion(1).cell.Ds * fphys ));
            
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ fluxP ] = pfieldP .* edge.nx;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            % The numerical at inner edge is given as $\hat Q = \left {\nabla_h p^- + \nabla_h p+\right }-\tau [p]$
            [ fluxS ] = edge.nx .* (LocalVariableM + LocalVariableP)./2 - edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.nx - InnerEdgefp .* edge.nx);
%             [ fluxS ] = edge.nx .* (LocalVariableM + LocalVariableP)./2 - edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.ny - InnerEdgefp .* edge.ny);
%             [ fluxS ] = edge.nx .* (pfieldM + pfieldP)./2 - edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.nx - InnerEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ fluxS ] = zeros(size(fluxM));
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            % The numerical at Newmann edge is given as $\hat Q = p_N$
            ind = (edge.ftype == enumBoundaryCondition.Newmann );
            fluxS(:,ind) = edge.nx(:,ind) .* edge.nx(:,ind) .* physClass.GradExt(:, ind);
            % The numerical at Dirichlet edge is given as $\hat Q = \nabla_h p - \tau (u^--u_D)\boldsymbol{n}$
            ind = (edge.ftype == enumBoundaryCondition.Dirichlet );
            fluxS(:,ind) = edge.nx(:,ind) .* (LocalVariableM(:,ind) - ...
                edge.nx(:,ind) .* obj.BoundaryEdgeTau(:,ind)./Prantl.*(BoundaryEdgefm(:,ind) - BoundaryEdgefp(:,ind)));            
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
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
            frhs = physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * pfield ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * pfield );
            
            LocalVariable = Kappa .* ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * fphys ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * fphys ) );
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ fluxP ] = pfieldP .* edge.ny;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.ny .* (LocalVariableM + LocalVariableP)./2 - ...
                edge.ny .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.ny - InnerEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ fluxS ] = zeros(size(fluxM));
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            ind = (edge.ftype == enumBoundaryCondition.Newmann );
            fluxS(:,ind) = edge.ny(:,ind) .*  edge.ny(:,ind) .* physClass.GradExt(:, ind);
            % The numerical at Dirichlet edge is given as $\hat Q = \nabla_h p - \tau (u^--u_D)\boldsymbol{n}$
            ind = (edge.ftype == enumBoundaryCondition.Dirichlet );
            fluxS(:,ind) = edge.ny(:,ind) .* (LocalVariableM(:,ind) - ...
                edge.ny(:,ind) .* obj.BoundaryEdgeTau(:,ind)./Prantl.*(BoundaryEdgefm(:,ind) - BoundaryEdgefp(:,ind)));            
%             [ fluxS ] = edge.nx .* LocalVariableM - edge.nx .* obj.BoundaryEdgeTau./Prantl .* ( BoundaryEdgefm .* edge.nx - BoundaryEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM,  fluxS);
        end
        
        function frhs = matCalculateMixedPartDerivTermX(obj, physClass, pfield, Kappa, fphys, Prantl, InnerEdgefm, InnerEdgefp, BoundaryEdgefm, BoundaryEdgefp)
            %> @brief Calculating the mixed derivative in x direction
            %> @detail this function is used to calculate the mixed
            %> derivative in x direction, i.e. $\frac{\partial \sigma}{\partial x}$,
            %> in the above formulation, $\sigma$ is $\Kappa\frac{\partial v}{\partial y}$
            %> @param[in] physClass The physical solver establised
            %> @param[in] pfield Value of the first order derivative of v,'$\frac{\partial v}{\partial y}$'
            %> this value has been contained obj.py, we just need to fetch it as we need
            %> @param[in] Kappa The diffusion coefficient, this is needed to
            %> calculate the numerical flux
            %> @param[in] fphys The physical field to be calculated, i.e. v in the presentation
            %> @param[in] Prantl The Prantl number, for conventional convection-diffusion
            %> problem, this value equals to 1, however, for salt and tempreture contained
            %> in SWE3d, this value equals to 1.1, this value is needed here, since it
            %> influence Kappa and the penalty parameter.
            %> @param[in] InnerEdgefm The local value of fphys at InnerEdge, used to
            %> calculate numerical flux
            %> @param[in] InnerEdgefp The adjacent value of fphys at InnerEdge, used
            %> to calculate numerical flux
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
            frhs = physClass.meshUnion(1).rx .* ( physClass.meshUnion(1).cell.Dr * pfield ) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds * pfield );
            
            LocalVariable = Kappa .* ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * fphys ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * fphys ) );
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ fluxP ] = pfieldP .* edge.nx;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.nx .* (LocalVariableM + LocalVariableP)./2 - ...
                edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.ny - InnerEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.nx .* LocalVariableM - edge.nx .* obj.BoundaryEdgeTau./Prantl .* ( BoundaryEdgefm .* edge.nx - BoundaryEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS);
        end
        
        function frhs = matCalculateMixedPartDerivTermY(obj, physClass, pfield, Kappa, fphys, Prantl, InnerEdgefm, InnerEdgefp, BoundaryEdgefm, BoundaryEdgefp)
            %> @brief Calculating the mixed derivative in y direction
            %> @detail this function is used to calculate the mixed
            %> derivative in y direction, i.e. $\frac{\partial \sigma}{\partial y}$,
            %> in the above formulation, $\sigma$ is $\Kappa\frac{\partial u}{\partial x}$
            %> @param[in] physClass The physical solver establised
            %> @param[in] pfield Value of the first order derivative of u,'$\frac{\partial u}{\partial x}$'
            %> this value has been contained obj.px, we just need to fetch it as we need
            %> @param[in] Kappa The diffusion coefficient, this is needed to
            %> calculate the numerical flux
            %> @param[in] fphys The physical field to be calculated, i.e. u in the presentation
            %> @param[in] Prantl The Prantl number, for conventional convection-diffusion
            %> problem, this value equals to 1, however, for salt and tempreture contained
            %> in SWE3d, this value equals to 1.1, this value is needed here, since it
            %> influence Kappa and the penalty parameter.
            %> @param[in] InnerEdgefm The local value of fphys at InnerEdge, used to
            %> calculate numerical flux
            %> @param[in] InnerEdgefp The adjacent value of fphys at InnerEdge, used
            %> to calculate numerical flux
            %> @param[in] BoundaryEdgefm The local value of fphys at BoundaryEdge
            %> @param[in] BoundaryEdgefp The adjacent value of fphys at BoundaryEdge
            frhs = physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * pfield ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * pfield );
            
            LocalVariable = Kappa .* ( physClass.meshUnion(1).rx .* ( physClass.meshUnion(1).cell.Dr * fphys ) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds * fphys ) );
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ fluxP ] = pfieldP .* edge.ny;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.ny .* (LocalVariableM + LocalVariableP)./2 - ...
                edge.ny .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.nx - InnerEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.ny .* LocalVariableM  - edge.ny .* obj.BoundaryEdgeTau./Prantl .* ( BoundaryEdgefm .* edge.ny - BoundaryEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS);
        end
        
        function matUpdatePenaltyParameter( obj, physClass, DiffusionCoefficient )
            %> this penalty parameter is calculated as $\tau=\frac{(D_p+1)(D_p+d)}{d}\frac{n_0}{2}\frac{A}{V}\miu$
            [ HnvM, HnvP ] = obj.matEvaluateSurfValue(physClass.meshUnion(1).InnerEdge, DiffusionCoefficient );
            
%             InnerEdgeA_fm = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(1,:)))', 1, physClass.meshUnion(1).Nz );
%             InnerEdgeA_fp = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(2,:)))', 1, physClass.meshUnion(1).Nz );
            InnerEdgeA_fm = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(1,:)))', 1, 1 );
            InnerEdgeA_fp = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(2,:)))', 1, 1 );
            InnerEdgeTau_fm = bsxfun(@times,  ( InnerEdgeA_fm(:) )',...
                ( physClass.meshUnion(1).cell.N + 1 )*(physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )/double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* HnvM);
            InnerEdgeTau_fp = bsxfun(@times,  ( InnerEdgeA_fp(:) )',...
                ( physClass.meshUnion(1).cell.N + 1 )*(physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )/double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* HnvP);
            obj.InnerEdgeTau = 1 * max( InnerEdgeTau_fm, InnerEdgeTau_fp );
            
%             BoundaryEdgeA_fm = repmat ( (physClass.mesh2d(1).BoundaryEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).BoundaryEdge.FToE(1,:)))', 1, physClass.meshUnion(1).Nz );
            BoundaryEdgeA_fm = repmat ( (physClass.mesh2d(1).BoundaryEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).BoundaryEdge.FToE(1,:)))', 1, 1 );
            [ Hnv, ~ ] = obj.matEvaluateSurfValue(physClass.meshUnion(1).BoundaryEdge, DiffusionCoefficient);
            obj.BoundaryEdgeTau = 1 * bsxfun(@times, ( BoundaryEdgeA_fm(:) )', ...
                ( physClass.meshUnion(1).cell.N + 1 )*( physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )./double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* Hnv);
        end
    end
end

