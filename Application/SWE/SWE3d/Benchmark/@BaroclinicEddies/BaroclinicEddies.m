classdef BaroclinicEddies < SWEBaroclinic3d
    %LOCKEXCHANGECASE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ChLength = 500000
        
        ChWidth = 160000
        
        finalTime = 201 * 86400
        
        H0 = 1000
        
        RPEfid
        
        RPE0
        
    end
    
    properties( Constant )
        hcrit = 0.01
    end
    
    methods
        %> For this case, the parameter is set following (Thetis, 2017), and we take Mx = 125, My = 40 and Mz = 40
        function obj = BaroclinicEddies(N, Nz, M, Mz)
            %LOCKEXCHANGECASE 构造此类的实例
            %   此处显示详细说明
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            obj.Cf{1} = 0.01*ones(size(obj.mesh2d.x));
                        
            obj.RPEfid = fopen('Result\LockExchangeCase\3d\RPEData.dat','w');
            
            obj.fphys{1}(:,:,13) = obj.matCalculateDensityField( obj.fphys{1} );
            
            obj.RPE0 = obj.matCalculateRPE( obj.fphys );
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
                fphys{m}(:,:,15) = 35*obj.H0;
                
                TBot = 10;
                TSurf = 20;
                TempData = obj.H0*(TBot + (TSurf-TBot)*((-1) - mesh3d.z)./(-1));
%                 fphys{m}(:,:,14) = obj.H0*(TBot + (TSurf-TBot)*((-1) - mesh3d.z)./(-1));
                Lx = 160000;
                y0 = 250000;
                k = 3;
                yA = 40000;
                deltaT = 1.2;
                deltaY = 40000;
                
                Index = mesh3d.y <= y0 - yA*sin(2*pi*k*mesh3d.x./Lx);
                TempData(Index) = TempData(Index) - deltaT;
                
                Index = mesh3d.y > y0 - yA*sin(2*pi*k*mesh3d.x./Lx) & mesh3d.y < y0 - yA*sin(2*pi*k*mesh3d.x./Lx) + deltaY;
                TempData(Index) = TempData(Index) - deltaT * ( 1 - (mesh3d.y(Index) - (y0 - yA*sin(2*pi*k*mesh3d.x(Index)./Lx) ))/40000);
                
                x2 = 110000;
                x3 = 130000;
                deltaT = 0.3;
                yw = y0 - yA/2*sin(pi*(mesh3d.x - x2)./(x3 - x2));
                Index = ( mesh3d.x >= x2 & mesh3d.x <= x3 & mesh3d.y >= yw - deltaY/2 & mesh3d.y <= yw + deltaY/2 );
                TempData(Index) = TempData(Index) + obj.H0 * deltaT*(1 - (mesh3d.y(Index) - yw(Index))/deltaY/2);
                fphys{m}(:,:,14) = TempData;
                
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0;
            end
        end
        
        function matUpdateOutputResult( obj, time, fphys2d, fphys )
            
            matUpdateOutputResult@SWEAbstract3d( obj, time, fphys2d, fphys );
                        
            RPE = obj.matCalculateRPE( fphys );
            
            fprintf(obj.RPEfid,'%12.8f  %12.8f\n', time, (RPE - obj.RPE0)/obj.RPE0 ); 
        end
        
        function RPE = matCalculateRPE( obj, fphys )
            
            aveRHO = obj.meshUnion(1).GetMeshAverageValue(...
                    fphys{1}(:,:,13) );
            [ sortedAveRHO, I ] = sort(aveRHO,'descend');
            %The mesh is Uniform
            aveZ = obj.meshUnion(1).GetMeshAverageValue(...
                    obj.meshUnion(1).z );
            sortedAveZ = sort(aveZ, 'ascend');
            sortedAveZ = sortedAveZ * obj.H0 + obj.H0 ;
            sortedLAV = obj.meshUnion(1).LAV(I);
            
            RPE = sum( obj.gra * sortedLAV .* sortedAveZ .* sortedAveRHO );
            
        end
        
        function matUpdateFinalResult( obj, ~, ~, ~ )
            
            fclose(obj.RPEfid);
            
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
           %doing nothing
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
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
%             option('GOTMSetupFile') = obj.GotmFile;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.0001;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.Constant;
            option('ConstantHorizontalEddyViscosityValue') = 1.0;
%             option('PhysicalSurfaceRoughnessLength') = 0.02;
%             option('PhysicalBottomRoughnessLength') = 0.0015;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.Beta;
            option('f0 for beta coriolis solver') = 1.2*10^(-4);
            option('beta for beta coriolis solver') = 0.0;
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
    [ 0, obj.ChLength ], [ 0, obj.ChWidth ], 125, 40, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
end 

