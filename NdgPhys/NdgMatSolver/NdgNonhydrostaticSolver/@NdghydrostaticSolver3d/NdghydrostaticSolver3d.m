classdef NdghydrostaticSolver3d
    %NDGHYDROSTATICSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        
        function obj = NdghydrostaticSolver3d(PhysClass)
            %Doing Nothing
        end
        
        function fphys = NdgConservativeNonhydrostaticUpdata(obj, physClass, fphys, fphys2d, deltatime)
            %Doing Nothing
        end
        
        function matGetBottomOldMomentum( obj, physClass, fphys )
            %Doing Nothing
        end
        
        function matClearGlobalMemory( obj )
            %Doing Nothing
        end
    end
    
end

