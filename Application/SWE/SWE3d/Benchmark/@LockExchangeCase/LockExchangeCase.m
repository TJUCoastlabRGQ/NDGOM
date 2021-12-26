classdef LockExchangeCase < SWEBaroclinic3d
    %LOCKEXCHANGECASE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ChLength = 64000
        
        ChWidth = 500
        
        finalTime = 61200
        
        H0 = 20
        
        MAXfid
        
        MINfid
        
        RPEfid
        
        RPE0
        
        Front95To96fid
        
        Front96To97fid
        
        Front97To98fid
        
        Front98To99fid
        
        Front99To00fid
    end
    
    properties( Constant )
        hcrit = 0.01
    end
    
    methods
        %> For this case, the parameter is set following (Thetis, 2017), and we take M = 128, Mz = 20
        function obj = LockExchangeCase(N, Nz, M, Mz)
            %LOCKEXCHANGECASE 构造此类的实例
            %   此处显示详细说明
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            obj.Cf{1} = 0*ones(size(obj.mesh2d.x));
            
            obj.MAXfid = fopen('Result\LockExchangeCase\3d\MaxData.dat','w');
            
            obj.MINfid = fopen('Result\LockExchangeCase\3d\MinData.dat','w');
            
            obj.RPEfid = fopen('Result\LockExchangeCase\3d\RPEData.dat','w');
            
            obj.Front95To96fid = fopen('Result\LockExchangeCase\3d\FrontData95To96.dat','w');
            
            obj.Front96To97fid = fopen('Result\LockExchangeCase\3d\FrontData96To97.dat','w');
            
            obj.Front97To98fid = fopen('Result\LockExchangeCase\3d\FrontData97To98.dat','w');
            
            obj.Front98To99fid = fopen('Result\LockExchangeCase\3d\FrontData98To99.dat','w');
            
            obj.Front99To00fid = fopen('Result\LockExchangeCase\3d\FrontData99To00.dat','w');
            
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
                Index = all(mesh3d.x<=obj.ChLength/2);
                fphys{m}(:,Index,15) = 35*obj.H0;
                Index = all(mesh3d.x>=obj.ChLength/2);
                fphys{m}(:,Index,15) = 35*obj.H0;
                Index = all(mesh3d.x<=obj.ChLength/2);
                fphys{m}(:,Index,14) = 5*obj.H0;
                Index = all(mesh3d.x>=obj.ChLength/2);
                fphys{m}(:,Index,14) = 30*obj.H0;
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0;
            end
        end
        
        function matUpdateOutputResult( obj, time, ~, fphys )
            
            fprintf(obj.MAXfid,'%12.8f  %12.8f\n', time, max(max(fphys{1}(:,:,13))) );
            
            fprintf(obj.MINfid,'%12.8f  %12.8f\n', time, min(min(fphys{1}(:,:,13))) );
            
            Index = (fphys{1}(:,:,13)>995 & fphys{1}(:,:,13)<=996);
            
            fprintf(obj.Front95To96fid,'%12.8f  %12.8f\n', time, max(max(obj.meshUnion.x(Index))) );
            
            Index = (fphys{1}(:,:,13)>996 & fphys{1}(:,:,13)<=997);
            
            fprintf(obj.Front96To97fid,'%12.8f  %12.8f\n', time, max(max(obj.meshUnion.x(Index))) ); 
            
            Index = (fphys{1}(:,:,13)>997 & fphys{1}(:,:,13)<=998);
            
            fprintf(obj.Front97To98fid,'%12.8f  %12.8f\n', time, max(max(obj.meshUnion.x(Index))) );
            
            Index = (fphys{1}(:,:,13)>998 & fphys{1}(:,:,13)<=999);
            
            fprintf(obj.Front98To99fid,'%12.8f  %12.8f\n', time, max(max(obj.meshUnion.x(Index))) );
            
            Index = (fphys{1}(:,:,13)>999 & fphys{1}(:,:,13)<=1000);
            
            fprintf(obj.Front99To00fid,'%12.8f  %12.8f\n', time, max(max(obj.meshUnion.x(Index))) ); 
            
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
            sortedAveZ = aveZ(I) * obj.H0 + obj.H0 ;
            sortedLAV = obj.meshUnion(1).LAV(I);
            
            RPE = sum( obj.gra * sortedLAV .* sortedAveZ .* sortedAveRHO );
            
        end
        
        function matUpdateFinalResult( obj, ~, ~, ~ )
            fclose(obj.MAXfid);
            
            fclose(obj.MINfid);
            
            fclose(obj.RPEfid);
            
            fclose(obj.Front95To96fid);
                    
            fclose(obj.Front96To97fid);
            
            fclose(obj.Front97To98fid);
            
            fclose(obj.Front98To99fid);
            
            fclose(obj.Front99To00fid);
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