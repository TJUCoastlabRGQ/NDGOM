classdef NdgSWEHorizDiffSolver < NdgHorizDiffSolver
    %NDGSWEHORIZDIFFSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function matEvaluateDiffRHS(obj, physClass, fphys)
            
            obj.matUpdateViscosity(physClass, fphys(:,:,1), fphys(:,:,2), fphys(:,:,4));
            obj.matUpdatePenaltyParameter( physClass, fphys(:,:,4) );
            %> $\Kappa = \nv * H$
            Kappa = obj.nv .* fphys(:,:,4);
            for i = 1:2
                obj.matCalculateAuxialaryVariable( physClass, fphys(:,:,physClass.varFieldIndex(i))./fphys(:,:,4), Kappa, i, ...
                    physClass.InnerEdgefm3d{1}(:,:,physClass.varFieldIndex(i))./physClass.InnerEdgefm3d{1}(:,:,4),...
                    physClass.InnerEdgefp3d{1}(:,:,physClass.varFieldIndex(i))./physClass.InnerEdgefp3d{1}(:,:,4), ...
                    physClass.BoundaryEdgefm3d{1}(:,:,physClass.varFieldIndex(i))./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                    physClass.BoundaryEdgefp3d{1}(:,:,physClass.varFieldIndex(i))./physClass.BoundaryEdgefp3d{1}(:,:,4));
            end
            for i = 3:physClass.Nvar
                obj.matCalculateAuxialaryVariable( physClass, fphys(:,:,physClass.varFieldIndex(i))./fphys(:,:,4), Kappa./obj.Prantl, i, ...
                    physClass.InnerEdgefm3d{1}(:,:,physClass.varFieldIndex(i))./physClass.InnerEdgefm3d{1}(:,:,4),...
                    physClass.InnerEdgefp3d{1}(:,:,physClass.varFieldIndex(i))./physClass.InnerEdgefp3d{1}(:,:,4), ...
                    physClass.BoundaryEdgefm3d{1}(:,:,physClass.varFieldIndex(i))./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                    physClass.BoundaryEdgefp3d{1}(:,:,physClass.varFieldIndex(i))./physClass.BoundaryEdgefp3d{1}(:,:,4));
            end
            %> this part is used to calculate $\frac{\partial}{\partial x}(2\nv_h H\frac{\partial u}{\partial x})
            %> + \frac{\partial}{\partial y}(\nv H(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}))$
            physClass.frhs{1}(:,:,1) = physClass.frhs{1}(:,:,1) + 2 * obj.matCalculatePartDerivTermX( physClass, obj.px(:,:,1),...
                Kappa, fphys(:,:,1)./fphys(:,:,4), 1, physClass.InnerEdgefm3d{1}(:,:,1)./physClass.InnerEdgefm3d{1}(:,:,4), ...
                physClass.InnerEdgefp3d{1}(:,:,1)./physClass.InnerEdgefp3d{1}(:,:,4), physClass.BoundaryEdgefm3d{1}(:,:,1)./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                physClass.BoundaryEdgefp3d{1}(:,:,1)./physClass.BoundaryEdgefp3d{1}(:,:,4) ) + ...
                obj.matCalculatePartDerivTermY( physClass, obj.py(:,:,1),...
                Kappa, fphys(:,:,1)./fphys(:,:,4), 1, physClass.InnerEdgefm3d{1}(:,:,1)./physClass.InnerEdgefm3d{1}(:,:,4), ...
                physClass.InnerEdgefp3d{1}(:,:,1)./physClass.InnerEdgefp3d{1}(:,:,4), physClass.BoundaryEdgefm3d{1}(:,:,1)./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                physClass.BoundaryEdgefp3d{1}(:,:,1)./physClass.BoundaryEdgefp3d{1}(:,:,4) ) + ...
                obj.matCalculatePartDerivTermY( physClass, obj.px(:,:,2),...
                Kappa, fphys(:,:,2)./fphys(:,:,4), 1, physClass.InnerEdgefm3d{1}(:,:,2)./physClass.InnerEdgefm3d{1}(:,:,4), ...
                physClass.InnerEdgefp3d{1}(:,:,2)./physClass.InnerEdgefp3d{1}(:,:,4), physClass.BoundaryEdgefm3d{1}(:,:,2)./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                physClass.BoundaryEdgefp3d{1}(:,:,2)./physClass.BoundaryEdgefp3d{1}(:,:,4) );
            
            %> this part is used to calculate $\frac{\partial}{\partial x}(\nv H(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}))+
            %> \frac{\partial}{\partial y}(2\nv_h H\frac{\partial v}{\partial y} $
            physClass.frhs{1}(:,:,2) = physClass.frhs{1}(:,:,2) +  obj.matCalculatePartDerivTermX( physClass, obj.py(:,:,1),...
                Kappa, fphys(:,:,1)./fphys(:,:,4), 1, physClass.InnerEdgefm3d{1}(:,:,1)./physClass.InnerEdgefm3d{1}(:,:,4), ...
                physClass.InnerEdgefp3d{1}(:,:,1)./physClass.InnerEdgefp3d{1}(:,:,4), physClass.BoundaryEdgefm3d{1}(:,:,1)./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                physClass.BoundaryEdgefp3d{1}(:,:,1)./physClass.BoundaryEdgefp3d{1}(:,:,4) ) + ...
                obj.matCalculatePartDerivTermX( physClass, obj.px(:,:,2),...
                Kappa, fphys(:,:,2)./fphys(:,:,4), 1, physClass.InnerEdgefm3d{1}(:,:,2)./physClass.InnerEdgefm3d{1}(:,:,4), ...
                physClass.InnerEdgefp3d{1}(:,:,2)./physClass.InnerEdgefp3d{1}(:,:,4), physClass.BoundaryEdgefm3d{1}(:,:,2)./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                physClass.BoundaryEdgefp3d{1}(:,:,2)./physClass.BoundaryEdgefp3d{1}(:,:,4) ) + ...
                2 * obj.matCalculatePartDerivTermY( physClass, obj.py(:,:,2),...
                Kappa, fphys(:,:,2)./fphys(:,:,4), 1, physClass.InnerEdgefm3d{1}(:,:,2)./physClass.InnerEdgefm3d{1}(:,:,4), ...
                physClass.InnerEdgefp3d{1}(:,:,2)./physClass.InnerEdgefp3d{1}(:,:,4), physClass.BoundaryEdgefm3d{1}(:,:,2)./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                physClass.BoundaryEdgefp3d{1}(:,:,2)./physClass.BoundaryEdgefp3d{1}(:,:,4) );
            
            for i = 3:physClass.Nvar
                physClass.frhs{1}(:,:,i) = physClass.frhs{1}(:,:,i) + obj.matCalculatePartDerivTermX( physClass, obj.px(:,:,i),...
                    Kappa./obj.Prantl, fphys(:,:,physClass.varFiledIndex(i))./fphys(:,:,4), obj.Prantl, ...
                    physClass.InnerEdgefm3d{1}(:,:,physClass.varFiledIndex(i))./physClass.InnerEdgefm3d{1}(:,:,4), ...
                    physClass.InnerEdgefp3d{1}(:,:,physClass.varFiledIndex(i))./physClass.InnerEdgefp3d{1}(:,:,4),...
                    physClass.BoundaryEdgefm3d{1}(:,:,physClass.varFiledIndex(i))./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                    physClass.BoundaryEdgefp3d{1}(:,:,physClass.varFiledIndex(i))./physClass.BoundaryEdgefp3d{1}(:,:,4) )+ ...
                    obj.matCalculatePartDerivTermY( physClass, obj.py(:,:,i),...
                    Kappa./obj.Prantl, fphys(:,:,physClass.varFiledIndex(i))./fphys(:,:,4), obj.Prantl, ...
                    physClass.InnerEdgefm3d{1}(:,:,physClass.varFiledIndex(i))./physClass.InnerEdgefm3d{1}(:,:,4), ...
                    physClass.InnerEdgefp3d{1}(:,:,physClass.varFiledIndex(i))./physClass.InnerEdgefp3d{1}(:,:,4),...
                    physClass.BoundaryEdgefm3d{1}(:,:,physClass.varFiledIndex(i))./physClass.BoundaryEdgefm3d{1}(:,:,4),...
                    physClass.BoundaryEdgefp3d{1}(:,:,physClass.varFiledIndex(i))./physClass.BoundaryEdgefp3d{1}(:,:,4) );
            end
        end
    end
    
    methods( Access = protected )
        function matUpdatePenaltyParameter( obj, physClass, Height )
            %> this penalty parameter is calculated as $\tau=\frac{(D_p+1)(D_p+d)}{d}\frac{n_0}{2}\frac{A}{V}\miu$
            [ HnvM, HnvP ] = obj.matEvaluateSurfValue(physClass.meshUnion(1).InnerEdge, obj.nv .* Height);
            
            InnerEdgeA_fm = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(1,:)))', 1, physClass.meshUnion(1).Nz );
            InnerEdgeA_fp = repmat( (physClass.mesh2d(1).InnerEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).InnerEdge.FToE(2,:)))', 1, physClass.meshUnion(1).Nz );
            InnerEdgeTau_fm = bsxfun(@times,  ( InnerEdgeA_fm(:) )',...
                ( physClass.meshUnion(1).cell.N + 1 )*(physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )/double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* HnvM);
            InnerEdgeTau_fp = bsxfun(@times,  ( InnerEdgeA_fp(:) )',...
                ( physClass.meshUnion(1).cell.N + 1 )*(physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )/double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* HnvP);
            obj.InnerEdgeTau = max( InnerEdgeTau_fm, InnerEdgeTau_fp );
            
            BoundaryEdgeA_fm = repmat( (physClass.mesh2d(1).BoundaryEdge.LAV./physClass.mesh2d(1).LAV(physClass.mesh2d(1).BoundaryEdge.FToE(1,:)))', 1, physClass.meshUnion(1).Nz );
            [ Hnv, ~ ] = obj.matEvaluateSurfValue(physClass.meshUnion(1).BoundaryEdge, obj.nv .* Height);
            obj.BoundaryEdgeTau = bsxfun(@times, ( BoundaryEdgeA_fm(:) )', ...
                ( physClass.meshUnion(1).cell.N + 1 )*( physClass.meshUnion(1).cell.N + ...
                double(physClass.meshUnion(1).type) )./double(physClass.meshUnion(1).type) * physClass.meshUnion(1).cell.Nface/2 .* Hnv);
        end
    end
    
end
