classdef NdgSlopeBottomWetDryTestWithWallBoundary < NdgFlatBottomWetDryTestWithWallBoundary
    %NDGSLOPEBOTTOMWETDRYTESTWITHWALLBOUNDARY 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgSlopeBottomWetDryTestWithWallBoundary(N, cellType)
            obj = obj@NdgFlatBottomWetDryTestWithWallBoundary(N, cellType);
        end
        
        function ExactGetFaceOuterValue = getExactGetFaceOuterValue(obj)
            %> In this case, the nonhydrostatic pressure is set to be one
            Nfp = obj.meshUnion(1).cell.Nfp;
            Nonhydro = obj.zGrad{1}(:,:,1);
            ExactGetFaceOuterValue = Nonhydro(obj.meshUnion(1).eidP);
            ExactGetFaceOuterValue(2*Nfp+1:3*Nfp,2) = -1;
            ExactGetFaceOuterValue(Nfp+1:2*Nfp,4) = -1;
            ExactGetFaceOuterValue(3*Nfp+1:4*Nfp,6) = -1;
            ExactGetFaceOuterValue(1:Nfp,8) = -1;
            ExactGetFaceOuterValue(:,5) = 0;
        end
    end
    
    methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell(obj.Nmesh, 1);
            fphys{1}(:,:,4) = obj.meshUnion(1).x;
            fphys{1}(:,:,1) = 4 - fphys{1}(:,:,4);
            fphys{1}(:,5,1) = 0;
        end
    end
    
end

