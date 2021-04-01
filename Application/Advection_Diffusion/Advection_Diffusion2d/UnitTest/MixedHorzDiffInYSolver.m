classdef MixedHorzDiffInYSolver < MixedHorzDiffSolver
    
    methods
        function obj = MixedHorzDiffInYSolver( physClass )
            obj = obj@MixedHorzDiffSolver( physClass );
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
                StiffMatrix = obj.matCalculateMixedPartDerivTermY( physClass, obj.px(:,:,i),...
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
                physClass.frhs{1}(:,:,i) = physClass.frhs{1}(:,:,i) + obj.matCalculateMixedPartDerivTermY( physClass, obj.px(:,:,i),...
                    Kappa, fphys{1}(:,:,physClass.varFieldIndex(i)), 1, IEfm(:,:,physClass.varFieldIndex(i)), ...
                    IEfp(:,:,physClass.varFieldIndex(i)), BEfm(:,:,physClass.varFieldIndex(i)),...
                    BEfp(:,:,physClass.varFieldIndex(i)));
            end
            
        end
    end
end
