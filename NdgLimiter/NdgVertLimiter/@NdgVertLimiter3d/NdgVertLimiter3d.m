classdef NdgVertLimiter3d < NdgAbstractLimiter
    %NDGVERTLIMITER3D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        
        function obj = NdgVertLimiter3d( mesh )
            obj = obj@NdgAbstractLimiter( mesh );
        end
        
        function fphys = matLimit( obj, fphys, fieldId )
            % calculate the cell averages
            Nmesh = obj.Nmesh;
            cvar = cell( Nmesh, 1 );
            for m = 1:Nmesh
                cvar{m} = obj.meshUnion(m).GetMeshAverageValue(...
                    fphys{m}(:,:,fieldId) );
                %The following way can be only used for constant Jacobi
%                 Data = obj.meshUnion.cell.V\fphys{m}(:,:,fieldId);
%                 data = Data(1,:) .* obj.meshUnion.cell.V(1) ;
%                 difference = data - cvar{m};
               fphys{m}(:,:,fieldId) = mxVertLimit3d( fphys{m}(:,:,fieldId), cvar{m},...
                   obj.meshUnion(m).EToE, inv(obj.meshUnion(m).cell.V1d),...
                   obj.meshUnion(m).cell.Npz, obj.meshUnion(m).cell.Nph, ...
                   obj.meshUnion(m).cell.V1d);
            end
        end
    end
    
end

