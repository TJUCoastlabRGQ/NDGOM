classdef MixedHorzDiffSolver < NdgHorizDiffSolver
    
    methods
        function obj = MixedHorzDiffSolver( physClass )
            obj = obj@NdgHorizDiffSolver( physClass );
        end
        
        function StiffMatrix = matEvaluateStiffMatrixInPointForm( obj, physClass, fphys)
            mesh = physClass.meshUnion(1);
            edge = mesh.InnerEdge;
            [ IEfm, IEfp ] = edge.matEvaluateSurfValue( fphys );
            
            edge = mesh.BoundaryEdge;
            [ BEfm, BEfp ] = edge.matEvaluateSurfValue( fphys );
            [ ~, BEfp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, BEfm, BEfp, physClass.fext{1} );
            Kappa = obj.nv;
            for i = 1:physClass.Nvar
                obj.matCalculateAuxialaryVariable( physClass, fphys{1}(:,:,physClass.varFieldIndex(i)), Kappa, i, ...
                    IEfm(:,:,physClass.varFieldIndex(i)),...
                    IEfp(:,:,physClass.varFieldIndex(i)), ...
                    BEfm(:,:,physClass.varFieldIndex(i)),...
                    BEfp(:,:,physClass.varFieldIndex(i)));
            end
            for i = 1:physClass.Nvar
%                 StiffMatrix = obj.matCalculateMixedPartDerivTermX( physClass, obj.py(:,:,i),...
%                     Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
%                     IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
%                     BEfp(:,:,physClass.varFieldIndex(i)));
                StiffMatrix = obj.matCalculatePartDerivTermX( physClass, obj.px(:,:,i),...
                    Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
                    IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
                    BEfp(:,:,physClass.varFieldIndex(i))) + ...
                    obj.matCalculatePartDerivTermY( physClass, obj.py(:,:,i),...
                    Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
                    IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
                    BEfp(:,:,physClass.varFieldIndex(i)));
            end
        end
        
        function matEvaluateDiffRHS(obj, physClass, fphys)
            mesh = physClass.meshUnion(1);
            edge = mesh.InnerEdge;
            [ IEfm, IEfp ] = edge.matEvaluateSurfValue( fphys );
            
            edge = mesh.BoundaryEdge;
            [ BEfm, BEfp ] = edge.matEvaluateSurfValue( fphys );
            
            [ ~, BEfp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, BEfm, BEfp, physClass.fext{1} );
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
%                 physClass.frhs{1}(:,:,i) = physClass.frhs{1}(:,:,i) + obj.matCalculateMixedPartDerivTermX( physClass, obj.py(:,:,i),...
%                     Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
%                     IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
%                     BEfp(:,:,physClass.varFieldIndex(i))) + obj.matCalculateMixedPartDerivTermY( physClass, obj.px(:,:,i),...
%                     Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
%                     IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
%                     BEfp(:,:,physClass.varFieldIndex(i)))+ obj.matCalculatePartDerivTermX( physClass, obj.px(:,:,i),...
%                     Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
%                     IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
%                     BEfp(:,:,physClass.varFieldIndex(i))) + ...
%                     obj.matCalculatePartDerivTermY( physClass, obj.py(:,:,i),...
%                     Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
%                     IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
%                     BEfp(:,:,physClass.varFieldIndex(i)));
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
%             [ fluxS ] = edge.nx .* ((LocalVariableM + LocalVariableP)./2 - ...
%                 edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm - InnerEdgefp ));
            [ fluxS ] = edge.nx .* ((pfieldM + pfieldP)./2 - ...
                edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm - InnerEdgefp ));

            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ fluxS ] = zeros(size(fluxM));
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            ind = (edge.ftype == enumBoundaryCondition.Newmann );
            fluxS(:,ind) = edge.nx(:,ind) .* edge.nx(:,ind) .* physClass.GradExt(:, ind);
            ind = (edge.ftype == enumBoundaryCondition.Dirichlet );
            fluxS(:,ind) = edge.nx(:,ind) .* (LocalVariableM(:,ind) - ...
                edge.nx(:,ind) .* obj.BoundaryEdgeTau(:,ind)./Prantl.*(BoundaryEdgefm(:,ind) - BoundaryEdgefp(:,ind)));
            fluxS(:,ind) = edge.nx(:,ind) .* (pfieldM(:,ind) - ...
                edge.nx(:,ind) .* obj.BoundaryEdgeTau(:,ind)./Prantl.*(BoundaryEdgefm(:,ind) - BoundaryEdgefp(:,ind)));
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM,  fluxS);
        end
        
        function frhs = matCalculateMixedPartDerivTermY(obj, physClass, pfield, Kappa, fphys, Prantl, InnerEdgefm, InnerEdgefp, BoundaryEdgefm, BoundaryEdgefp)
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
                edge.ny .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.ny - InnerEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, ~] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ fluxS ] = zeros(size(fluxM));
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            ind = (edge.ftype == enumBoundaryCondition.Newmann );
            fluxS(:,ind) = edge.ny(:,ind) .* edge.ny(:,ind) .* physClass.GradExt(:, ind);
            ind = (edge.ftype == enumBoundaryCondition.Dirichlet );
            fluxS(:,ind) = edge.ny(:,ind) .* (LocalVariableM(:,ind) - ...
                edge.ny(:,ind) .* obj.BoundaryEdgeTau(:,ind)./Prantl.*(BoundaryEdgefm(:,ind) - BoundaryEdgefp(:,ind)));
            frhs = frhs - edge.matEvaluateStrongFormEdgeRHS( fluxM,  fluxS);
        end
        
    end
    
end