classdef LockExchangeCase < SWEBaroclinic3d
    %LOCKEXCHANGECASE �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
        ChLength = 2
        
        ChWidth = 0.01
        
        finalTime = 15
        
        H0 = 0.3
    end
    
    properties( Constant )
        hcrit = 0.01
    end
    
    methods
        function obj = LockExchangeCase(N, Nz, M, Mz)
            %LOCKEXCHANGECASE ��������ʵ��
            %   �˴���ʾ��ϸ˵��
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            obj.Cf{1} = 0*ones(size(obj.mesh2d.x));
        end
       
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                Index = all(mesh3d.x<=obj.ChLength/2);
                fphys{m}(:,Index,15) = 5*obj.H0;
                Index = all(mesh3d.x>=obj.ChLength/2);
                fphys{m}(:,Index,15) = 0.2*obj.H0;
                fphys{m}(:,:,14) = 20*obj.H0;
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
           
           h3d = zeros(size(obj.fext3d{1}(:,:,1)));
           h2d = zeros(size(obj.fext2d{1}(:,:,1)));
           Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
           h3d(:,Index) = -obj.meshUnion(1).BoundaryEdge.xb(:,Index) *10^(-5)  + 15;
           obj.fext3d{1}(:,:,3) = h3d;
           
           Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
           h2d(:,Index) = -obj.meshUnion(1).mesh2d.BoundaryEdge.xb(:,Index) *10^(-5)  + 15;
           obj.fext2d{1}(:,:,3) = h2d;
        end        
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 7500;
            option('startTime') = 0.0;
            option('finalTime') = obj.finalTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = obj.finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            %             option('EddyViscosityType') = enumEddyViscosity.GOTM;
            %             option('GOTMSetupFile') = obj.GotmFile;
            %             option('equationType') = enumDiscreteEquation.Strong;
            %             option('integralType') = enumDiscreteIntegral.QuadratureFree;
            %             option('outputType') = enumOutputFile.VTK;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
%             option('GOTMSetupFile') = obj.GotmFile;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 0;
%             option('PhysicalSurfaceRoughnessLength') = 0.02;
%             option('PhysicalBottomRoughnessLength') = 0.0015;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
        end
        
    end    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [ 0, obj.ChWidth ], M, 1, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
end