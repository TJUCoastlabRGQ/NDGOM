classdef NdgVertDiffSolver < AbstractDiffSolver
    %NDGVERTDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        BoundaryEdgeType = 'Dirichlet'
    end
    
    methods
        
        function obj = NdgVertDiffSolver( physClass )
            obj = obj@AbstractDiffSolver( physClass );
            if physClass.option.isKey('AdvDiffVerticalDiffusionType')
                if physClass.option.isKey('AdvDiffConstantVerticalDiffusionValue')
                    value = physClass.getOption('AdvDiffConstantVerticalDiffusionValue');
                    obj.nv = value * ones(size(physClass.meshUnion(1).x));
                    fprintf('Value of the constant vertical diffusion coefficient is set to be: %f\n',value);
                    obj.matUpdatePenaltyParameter(  physClass, obj.nv );
                end
            end
            if physClass.option.isKey('BottomBoundaryEdgeType')
                obj.BoundaryEdgeType = char(physClass.getOption('BottomBoundaryEdgeType'));
            end
            fprintf('The bottom boundary condition for momentum is: %s\n',obj.BoundaryEdgeType);
            
            obj.matClearGlobalMemory( );
        end
        
        function matClearGlobalMemory( obj )
            clear mxUpdateImplicitRHS;
        end
        
        %> @brief Calculating the right hand side corresponding to the vertical diffusion term and
        %> return the physical field with vertical diffusion considered
        %> @detail this function is used to calculate the right hand side corresponding to the vertical
        %> diffusion term and return the updated physical field at each Runge-Kutta time stage
        %> @param[in] physClass The physical solver establised
        %> @param[in] Height The water depth
        %> @param[in] ImplicitParameter The implicit parameter at the corresponding IMEXRK stage
        %> @param[in] dt The time step
        %> @param[out] fphys The physical field with vertical diffusion
        %> considered
        %> Input parameter changed on 20211231 to consider the hu and hv
        %> field, since we need it when we treat the bottom boundary implicitly. 
        function  fphys = matUpdateImplicitVerticalDiffusion( obj, physClass, SystemRHS, ImplicitParameter, dt, intRK, Stage, huv3d, h2d )
            fphys = obj.matCalculateImplicitRHS( physClass, obj.nv, SystemRHS, ImplicitParameter, dt, intRK, Stage, huv3d, h2d);
        end
        
        function matUpdateViscosity(obj, ~, ~, ~, ~, ~, ~ )
            %doing nothing
        end        
    end
    
    methods( Access = protected )
        %> Input parameter changed on 20211231 to consider the hu and hv
        %> field, since we need it when we treat the bottom boundary implicitly.  
        fphys  = matCalculateImplicitRHS( obj, physClass, DiffusionCoefficient, Height, SystemRHS, ImplicitParameter, dt, intRK, Stage, huv3d, h2d);
        
        function matUpdatePenaltyParameter( obj, physClass, DiffusionCoefficient )
            %> @brief Evaluating the penalty parameter used to penalize the jump between adjacet cell used in IPDG for second order operator
            %>@detail In this version, the Interior Penalty Discontinuous Galerkin(IPDG) method is used to treat
            %> the second order diffusion operator. To do so, the penalty parameter is calculated according to
            %> [1] Shahbazi K. An explicit expression for the penalty parameter of the interior penalty method[J].
            %> Journal of Computational Physics, 2005, 205(2): 401-407.
            %> [2] Pestiaux A. Parameterization of subgrid-scale processes in finite element sea ice-ocean models[D].
            %> UCL-Université Catholique de Louvain, 2015. pg:28.
            %> The formula is '$\tau=\frac{(D_p+1)(D_p+d)}{d}\frac{n_0}{2}\frac{A}{V}\miu$'
            %> @param[in] physClass The physical solver establised
            %> @param[in] DiffusionCoefficient The diffusion coefficient
            BotEidM   = physClass.meshUnion(1).cell.Fmask(physClass.meshUnion(1).cell.Fmask(:,end-1)~=0,end-1);
            UpEidM     = physClass.meshUnion(1).cell.Fmask(physClass.meshUnion(1).cell.Fmask(:,end)~=0,end);
            %             obj.tau = zeros( physClass.meshUnion(1).Nz+1, physClass.mesh2d(1).K );
            obj.tau = zeros( numel(BotEidM), physClass.mesh2d(1).K * ( physClass.meshUnion(1).Nz+1 ) );
            P = physClass.mesh2d(1).cell.N;
            %> for tri-prisms, number of faces is 5, for quad-prism, number of face is 6
            n0 = physClass.meshUnion(1).cell.Nface;
            %> here Nz stands for ratio between area of surface and volume of the studied cell
            Nz = physClass.meshUnion(1).Nz;
            for i = 1:physClass.mesh2d(1).K
                %> The surface most face for each column
                %     Tau(1,i) = (P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(UpEidM, (i-1)*Nz+1));
                for j = 2:Nz
                    %                     obj.tau(j,i) = (P+1)*(P+3)/3*n0/2*Nz*max(max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1)),...
                    %                         max(DiffusionCoefficient(UpEidM, (i-1)*Nz+j)));
                    obj.tau(:, (i-1)*( physClass.meshUnion(1).Nz+1 ) + j ) = (P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1),...
                        DiffusionCoefficient(UpEidM, (i-1)*Nz+j));
                end
                %> The bottom most face for each column
                %                 obj.tau(Nz+1,i) = (P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz));
                obj.tau(:, (i-1)*( physClass.meshUnion(1).Nz+1 ) + Nz + 1 ) = (P+1)*(P+3)/3*n0/2*Nz*DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz);
            end
        end
        
    end
    
end

