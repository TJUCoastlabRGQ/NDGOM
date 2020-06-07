classdef NdgHorizDiffSolver < AbstractDiffSolver
    %NDGHORIZDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties( Access = protected )
        px
        py
        InnerEdgeTau
        BoundaryEdgeTau
    end
    
    methods( Access = protected )
        
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
            
            obj.px(:,:,VarIndex) = Kappa .* ( physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr *  fphys ) + ...
                physClass.meshUnion(1).sx .* ( physClass.meshUnion(1).cell.Ds *  fphys ));
            
            obj.py(:,:,VarIndex) = Kappa .* ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * fphys ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * fphys ));
            
            edge = physClass.meshUnion(1).InnerEdge;
            %             [fM, fP] = obj.matEvaluateSurfValue(edge, fphys );
            [KappaM, KappaP] = obj.matEvaluateSurfValue(edge, Kappa );
            obj.px(:,:,VarIndex) = obj.px(:,:,VarIndex) - obj.matEvaluateStrongFormEdgeRHS(edge, ...
                InnerEdgefm .* edge.nx, InnerEdgefp .* edge.nx , 0.5 * (InnerEdgefm + InnerEdgefp) .* edge.nx, KappaM, KappaP);
            obj.py(:,:,VarIndex) = obj.py(:,:,VarIndex) - obj.matEvaluateStrongFormEdgeRHS(edge, ...
                InnerEdgefm .* edge.ny, InnerEdgefp .* edge.ny , 0.5 * (InnerEdgefm + InnerEdgefp) .* edge.ny, KappaM, KappaP);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            %             [fM, fP] = obj.matEvaluateSurfValue(edge, fphys );
            [KappaM, KappaP] = obj.matEvaluateSurfValue(edge, Kappa );
            obj.px(:,:,VarIndex) = obj.px(:,:,VarIndex) - obj.matEvaluateStrongFormEdgeRHS(edge, ...
                BoundaryEdgefm .* edge.nx, BoundaryEdgefp .* edge.nx , 0.5 * (BoundaryEdgefm + BoundaryEdgefp) .* edge.nx, KappaM, KappaP);
            obj.py(:,:,VarIndex) = obj.py(:,:,VarIndex) - obj.matEvaluateStrongFormEdgeRHS(edge, ...
                BoundaryEdgefm .* edge.ny, BoundaryEdgefp .* edge.ny , 0.5 * (BoundaryEdgefm + BoundaryEdgefp) .* edge.ny, KappaM, KappaP);
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
            LocalVariable = Kappa .* ( physClass.meshUnion(1).rx .* (physClass.meshUnion(1).cell.Dr * fphys ) + ...
                physClass.meshUnion(1).sx .* (physClass.meshUnion(1).cell.Ds * fphys ));
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ fluxP ] = pfieldP .* edge.nx;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.nx .* (LocalVariableM + LocalVariableP)./2 - edge.nx .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.nx - InnerEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.nx;
            [ fluxP ] = pfieldP .* edge.nx;
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.nx .* LocalVariableM - edge.nx .* obj.BoundaryEdgeTau./Prantl .* ( BoundaryEdgefm .* edge.nx - BoundaryEdgefp .* edge.nx);
            frhs = frhs - edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS);
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
            LocalVariable = Kappa .* ( physClass.meshUnion(1).ry .* ( physClass.meshUnion(1).cell.Dr * fphys ) + ...
                physClass.meshUnion(1).sy .* ( physClass.meshUnion(1).cell.Ds * fphys ) );
            edge = physClass.meshUnion(1).InnerEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ fluxP ] = pfieldP .* edge.ny;
            [ LocalVariableM, LocalVariableP ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.ny .* (LocalVariableM + LocalVariableP)./2 - ...
                edge.ny .* obj.InnerEdgeTau./Prantl .* ( InnerEdgefm .* edge.ny - InnerEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS);
            
            edge = physClass.meshUnion(1).BoundaryEdge;
            [pfieldM, pfieldP] = obj.matEvaluateSurfValue(edge, pfield);
            [ fluxM ] = pfieldM .* edge.ny;
            [ fluxP ] = pfieldP .* edge.ny;
            [ LocalVariableM, ~ ] = obj.matEvaluateSurfValue(edge, LocalVariable );
            [ fluxS ] = edge.ny .* LocalVariableM  - edge.ny .* obj.BoundaryEdgeTau./Prantl .* ( BoundaryEdgefm .* edge.ny - BoundaryEdgefp .* edge.ny);
            frhs = frhs - edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS);
        end
    end
    
end

