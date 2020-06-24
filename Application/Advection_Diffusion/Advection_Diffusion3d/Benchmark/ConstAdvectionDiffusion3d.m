classdef ConstAdvectionDiffusion3d < Adv_DiffAbstract3d
    %CONSTADVECTIONDIFFUSION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        miu = 0.01
    end
    
    methods
        function obj = ConstAdvectionDiffusion3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );            
        end
    end

    methods ( Access = protected )
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh3d = obj.mesh3d(m);
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                fphys{m}(:,:,1) = sin(2*pi*0)*sin(2*pi*mesh3d.x).*sin(pi*mesh3d.y).*...
                    sin(2*pi*mesh3d.z);
                fphys{m}(:,:,2) = obj.u0 * ones(size(mesh3d.x));
                fphys{m}(:,:,3) = obj.v0 * ones(size(mesh3d.x));
                fphys{m}(:,:,4) = obj.w0 * ones(size(mesh3d.x));
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            obj.BoundaryEdgefp3d{1} = sin(2*pi*time)*sin(2*pi*obj.mesh3d.BoundaryEdge.x).*...
                sin(pi*obj.mesh3d.BoundaryEdge.y).*sin(2*pi*obj.mesh3d.BoundaryEdge.z);
            obj.SurfaceBoundaryEdgefp3d{1} = sin(2*pi*time)*sin(2*pi*obj.mesh3d.SurfaceBoundaryEdge.x).*...
                sin(pi*obj.mesh3d.SurfaceBoundaryEdge.y).*sin(2*pi*obj.mesh3d.SurfaceBoundaryEdge.z);    
            obj.BottomBoundaryEdgefp3d{1} = sin(2*pi*time)*sin(2*pi*obj.mesh3d.BottomBoundaryEdge.x).*...
                sin(pi*obj.mesh3d.BottomBoundaryEdge.y).*sin(2*pi*obj.mesh3d.BottomBoundaryEdge.z);
            obj.SurfBoundNewmannDate{1}(:,:,1) = 2 * pi * obj.miu * sin(2*pi*time)*sin(2*pi*obj.mesh3d.SurfaceBoundaryEdge.x).*...
                sin(pi*obj.mesh3d.SurfaceBoundaryEdge.y).*cos(2*pi*obj.mesh3d.SurfaceBoundaryEdge.z) .* 1;
            obj.BotBoundNewmannDate{1}(:,:,1) = 2 * pi * obj.miu * sin(2*pi*time)*sin(2*pi*obj.mesh3d.BottomBoundaryEdge.x).*...
                sin(pi*obj.mesh3d.BottomBoundaryEdge.y).*cos(2*pi*obj.mesh3d.BottomBoundaryEdge.z) .* (-1);            
%         BotBoundNewmannDate
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 50;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('AdvDiffVerticalEddyViscosityType') = enumVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumHorizontalEddyViscosity.Constant;
            option('AdvDiffHorizontalEddyViscosityValue') = 100;
        end
        
    end    
end

