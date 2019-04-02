classdef SWEMSBarotropic3d < SWEBarotropic3d
    %SWEMSBAROTROPIC3D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
        Solver2d
    end
    
    methods
        function  obj = SWEMSBarotropic3d(  )
            %             obj.Solver2d = SWEPreBlanaced2d();
            %             obj.Solver2d.initPhysFromOptions( obj.mesh2d );
            obj.Nfield2d = 5; %[H HU HV eta b]
            obj.Nvar2d = 3;
            obj.varFieldIndex2d = [1 2 3];
           
        end
    end
    
    methods  ( Access = protected )
        fphys = matEvaluateCorrectedMomentum( obj, mesh3d, fphys, fphys2d );
        matEvaluateRK45( obj );
        matEvaluateRHS( obj, fphys2d, fphys );
    end
end

