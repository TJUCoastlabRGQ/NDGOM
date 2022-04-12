classdef IdealEstuaryCirculation < SWEBaroclinic3d
    %IDEALESTUARYCIRCULATION 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ChLength = 100000
        
        ChWidth = 2000
        
        H0 = 50
        
        finalTime = 16*24*3600
        
        GotmFile = fullfile([pwd,'/Application/SWE/SWE3d/Benchmark/@IdealEstuaryCirculation/gotmturb.nml'])
        
    end
    
    properties ( Constant )
        hcrit = 0.001;
    end
    
    methods
        function obj = IdealEstuaryCirculation(N, Nz, M, Mz)
            
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            
            obj.Nfield = 18;
            
            obj.fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy','w', 'hw','hc','rho', 'hT', 'hS', 'Tke', 'Eps','nvh'};
                        
            obj.outputFieldOrder3d = [1 2 3 4 5 6 14 15 16 17 18];
            
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
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
                fphys2d{m}(:,:,1) = -5/obj.ChLength*mesh2d.x + 7.5 + sin(pi*mesh2d.x/obj.ChLength);
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                %> For water depth
                fphys{m}(:,:,4) = -5/obj.ChLength*mesh3d.x + 7.5 + sin(pi*mesh3d.x/obj.ChLength);
                %> For temperature
                fphys{m}(:,:,14) = 10*fphys{m}(:,:,4);
                %> For salinity
                Index = all(mesh3d.x<=-20000);
                fphys{m}(:,Index,15) = 30 * fphys{m}(:,Index,4);
                Index = all(mesh3d.x>=-20000 & mesh3d.x<=30000);
                fphys{m}(:,Index,15) = (-30/50000*mesh3d.x(:,Index) + 90/5) .* fphys{m}(:,Index,4);
%                 fphys2d{m}(:, :, 4) = -1 * fphys2d{m}(:,:,1);
                fphys2d{m}(:, :, 4) = 5/obj.ChLength*mesh2d.x - 7.5;
            end
        end
        
        function matUpdateExternalField( obj, time, ~, fphys )
            hu3d = zeros(size(obj.fext3d{1}(:,:,1)));
            hs3d = zeros(size(obj.fext3d{1}(:,:,1)));
            ht3d = zeros(size(obj.fext3d{1}(:,:,1)));
            hu2d = zeros(size(obj.fext2d{1}(:,:,1)));
            %> The seaward boundary be set to ClampedVel
            Index = ( (obj.meshUnion(1).BoundaryEdge.ftype)' == enumBoundaryCondition.ClampedVel & all(obj.meshUnion(1).BoundaryEdge.xb == -obj.ChLength/2));
%             Momentum = -0.08*10 + 0.4*5*sin(2*pi*time/(12*3600));
            Momentum = -0.08*5 + 0.4*10*sin(2*pi*time/(12*3600));
            hu3d(:,Index) = Momentum;
            [fm, ~] = obj.meshUnion(1).BoundaryEdge.matEvaluateSurfValue(fphys);
            if ( Momentum >= 0)
                hs3d(:,Index) = fm(:,Index,4).*30;
            else
                hs3d(:,Index) = fm(:,Index,15);
            end
            ht3d(:,Index) = fm(:,Index,4).*10;
            %> The riverward boundary be set to ClampedVel
            Index = ( (obj.meshUnion(1).BoundaryEdge.ftype)' == enumBoundaryCondition.ClampedVel & all(obj.meshUnion(1).BoundaryEdge.xb == obj.ChLength/2));
            hu3d(:,Index) = -0.08*5;
            hs3d(:,Index) = fm(:,Index,4).*0;
            ht3d(:,Index) = fm(:,Index,4).*10;
            
            obj.fext3d{1}(:,:,1) = hu3d;
            %> here, index of the variable are arranged accoring to the number of the variable, not the index of the field
            obj.fext3d{1}(:,:,4) = ht3d;
%             obj.fext3d{1}(:,:,4) = fm(:,:,14);
            obj.fext3d{1}(:,:,5) = hs3d;
            %> For the 2d part
            Index = ( (obj.meshUnion(1).mesh2d.BoundaryEdge.ftype)' == enumBoundaryCondition.ClampedVel & all(obj.meshUnion(1).mesh2d.BoundaryEdge.xb == -obj.ChLength/2));
            hu2d(:,Index) = Momentum;
            Index = ( (obj.meshUnion(1).mesh2d.BoundaryEdge.ftype)' == enumBoundaryCondition.ClampedVel & all(obj.meshUnion(1).mesh2d.BoundaryEdge.xb == obj.ChLength/2));
            hu2d(:,Index) = -0.08*5;
            obj.fext2d{1}(:,:,1) = hu2d;
            
        end
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 800;
            option('startTime') = 0.0;
            option('finalTime') = obj.finalTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = obj.finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('EosType') = enumEOSType.Jackett05;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.0001;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('PhysicalSurfaceRoughnessLength') = 0.003;
            option('PhysicalBottomRoughnessLength') = 0.005;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.None;
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

% bctype = [ ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.ClampedVel, ...
%     enumBoundaryCondition.ClampedVel ];
bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], 1*[ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), 1*ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
end

