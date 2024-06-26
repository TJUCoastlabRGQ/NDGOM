classdef NdgVertLimiter2d < NdgVertLimiter
    
    methods
        function obj = NdgVertLimiter2d( mesh )
            obj = obj@NdgVertLimiter( mesh );
        end% func
        
        function fphys = matLimit( obj, physClass, fphys, fieldId )
            [ fvert, fvmin, fvmax, cvar ] = matEvaluateVertAverage( obj, fphys, fieldId );
            [ fvmax, fvmin ] = physClass.matEvaluateBoundaryVertScope( fphys, fvmax, fvmin );
            for m = 1:obj.Nmesh
                for i=1:numel(fieldId)
                fphys{m}(:,:,fieldId(i)) = mxVertLimit2d( ...
                    fphys{m}(:,:,fieldId(i)), ...
                    obj.meshUnion(m).x, ...
                    obj.meshUnion(m).y, ...
                    obj.meshUnion(m).xc, ...
                    obj.meshUnion(m).yc, ...
                    obj.meshUnion(m).vx, ...
                    obj.meshUnion(m).vy, ...
                    fvert{m}(:,:,i), ...
                    fvmin{m}(:,:,i), ...
                    fvmax{m}(:,:,i), ...
                    cvar{m}(:,:,i), ...
                    obj.meshUnion(m).EToV, ...
                    obj.meshUnion(m).cell.Fmask ...
                    );
                end
            end
        end% func
    end% methods
    
end

