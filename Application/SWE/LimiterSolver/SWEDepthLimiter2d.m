classdef SWEDepthLimiter2d < handle
    
    methods
        
        function fphys = apply( obj, phys, fphys )
            fphys = phys.limiter.matLimit( phys, fphys, phys.varFieldIndex );
%             fphys = phys.limiter.matLimit( phys, fphys, 1 );
%             fphys = phys.limiter.matLimit( phys, fphys, 2 );
%             fphys = phys.limiter.matLimit( phys, fphys, 3 );
        end
    end% methods
    
end

