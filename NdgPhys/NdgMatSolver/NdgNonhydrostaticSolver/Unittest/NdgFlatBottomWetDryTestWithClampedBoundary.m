classdef NdgFlatBottomWetDryTestWithClampedBoundary < NdgFlatBottomWetDryTestWithWallBoundary
    %UNTITLED �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgFlatBottomWetDryTestWithClampedBoundary(N, cellType)
            obj = obj@NdgFlatBottomWetDryTestWithWallBoundary(N, cellType);
        end
        
        function ExactReverseEidBoundaryType = getExactReverseEidBoundaryType(obj)
            ExactReverseEidBoundaryType = ones(size(obj.meshUnion(1).eidM));
        end
        
        function ExactEidBoundaryType = getExactEidBoundaryType(obj)
            Nfp = obj.meshUnion(1).cell.Nfp;
            ExactEidBoundaryType = ones(size(obj.meshUnion(1).eidM));
            ExactEidBoundaryType(1:Nfp,1) = -1;
            ExactEidBoundaryType(3*Nfp + 1 :4*Nfp,1) = -1;
            ExactEidBoundaryType(1:Nfp,2) = -1;
            ExactEidBoundaryType(1:2*Nfp,3) = -1;
            ExactEidBoundaryType(3*Nfp + 1 :4*Nfp,4) = -1;
            ExactEidBoundaryType(Nfp + 1 :2*Nfp,6) = -1;
            ExactEidBoundaryType(2*Nfp + 1 :4*Nfp,7) = -1;
            ExactEidBoundaryType(2*Nfp + 1 :3*Nfp,8) = -1;
            ExactEidBoundaryType(Nfp + 1 :3*Nfp,9) = -1;
        end
        
        function ExactGetFaceOuterValue = getExactGetFaceOuterValue(obj)
            Nfp = obj.meshUnion(1).cell.Nfp;
            Nonhydro = obj.fphys{1}(:,:,1);
            ExactGetFaceOuterValue = Nonhydro(obj.meshUnion(1).eidP);
            ExactGetFaceOuterValue(1:Nfp,1) = -1;
            ExactGetFaceOuterValue(3*Nfp+1:4*Nfp,1) = -1;
            ExactGetFaceOuterValue(1:Nfp,2) = -1;
            ExactGetFaceOuterValue(1:2*Nfp,3) = -1;
            ExactGetFaceOuterValue(3*Nfp+1:4*Nfp,4) = -1;
            ExactGetFaceOuterValue(Nfp+1:2*Nfp,6) = -1;
            ExactGetFaceOuterValue(2*Nfp+1:4*Nfp,7) = -1;
            ExactGetFaceOuterValue(2*Nfp+1:3*Nfp,8) = -1;
            ExactGetFaceOuterValue(Nfp+1:3*Nfp,9) = -1;
            ExactGetFaceOuterValue(2*Nfp+1:3*Nfp,2) = -1;
            ExactGetFaceOuterValue(Nfp+1:2*Nfp,4) = -1;
            ExactGetFaceOuterValue(3*Nfp+1:4*Nfp,6) = -1;
            ExactGetFaceOuterValue(1:Nfp,8) = -1;
            ExactGetFaceOuterValue(:,5) = 0;
        end  
        
        function ExactFluxterm = GetExactFluxterm(obj, mesh)
            Nfp = mesh.cell.Nfp;
            ExactFluxterm = ones(size(mesh.eidM));
            ExactFluxterm(2*Nfp+1:3*Nfp,2) = 0;
            ExactFluxterm(Nfp+1:2*Nfp,4) = 0;
            ExactFluxterm(3*Nfp+1:4*Nfp,6) = 0;
            ExactFluxterm(1:Nfp,8) = 0;    
            ExactFluxterm(:,5) = 0; 
        end
    end
    
    methods(Access = protected)
        function [ mesh ] = makeUniformMesh(obj, N, type)
            bctype = [...
                NdgEdgeType.Clamped, ...
                NdgEdgeType.Clamped, ...
                NdgEdgeType.Clamped, ...
                NdgEdgeType.Clamped];
            
            if (type == NdgCellType.Tri)
                mesh = makeUniformTriMesh(N, [0, 3], [0, 3], 3, 3, bctype);
            elseif(type == NdgCellType.Quad)
                mesh = makeUniformQuadMesh(N, [0, 3], [0, 3], 3, 3, bctype);
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func        
    end
    
end

