classdef WindDrivenFlow < SWEBarotropic3d
    %WINDDRIVENFLOW 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 2000;
        ChWidth = 300;
        %> channel depth
        H0 = 10;
        %> start time
        startTime = 0;
        %> final time
        finalTime = 3000;
        hcrit = 0.001;
    end
    
    properties
        dt
    end    
    
    methods
        function obj = WindDrivenFlow( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder = [ 1 2 3 4];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
%             obj.dt = 0.02;
            obj.Cf{1} = 0.0025*ones(size(obj.mesh2d(1).x));
            obj.WindTaux{1} = 1.5/1000*ones(size(obj.mesh2d(1).x));
            obj.WindTauy{1} = zeros(size(obj.mesh2d(1).y));            
%             obj.Cf{1} = 0.0025/1000;
        end
        
        matExplicitRK222(obj);
        
        
        EntropyAndEnergyCalculation(obj);
        
        AnalysisResult2d( obj );
        AnalysisResult3d( obj );
        
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
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;                
                %water depth
                fphys2d{m}(:,:,1) = obj.H0*ones(mesh2d.cell.Np, mesh2d.K);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 21600;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;                  
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('EddyViscosityType') = enumEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantEddyViscosityValue') = 0.01;
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
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition( obj, mesh2d, mesh3d );
end


function [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition( obj, mesh2d, mesh3d )
%we note that FTOF property is not corrected here, and this may be needed
%for non-hydrostatic nodel
%%periodic boundary condition for this case
LOrder2d = find(all(mesh2d.BoundaryEdge.xb == -obj.ChLength/2));
LOrder3d = find(all(mesh3d.BoundaryEdge.xb == -obj.ChLength/2));
ROrder2d = find(all(mesh2d.BoundaryEdge.xb == obj.ChLength/2));
ROrder3d = find(all(mesh3d.BoundaryEdge.xb == obj.ChLength/2));

%Ne part
mesh2d.BoundaryEdge.Ne = mesh2d.BoundaryEdge.Ne - numel(LOrder2d) - numel(ROrder2d);
mesh3d.BoundaryEdge.Ne = mesh3d.BoundaryEdge.Ne - numel(LOrder3d) - numel(ROrder3d);
mesh2d.InnerEdge.Ne = mesh2d.InnerEdge.Ne + numel(LOrder2d);
mesh3d.InnerEdge.Ne = mesh3d.InnerEdge.Ne + numel(LOrder3d);

%ftype part
mesh2d.BoundaryEdge.ftype([LOrder2d,ROrder2d])=[];
mesh3d.BoundaryEdge.ftype([LOrder3d,ROrder3d])=[];

%FToE part
InnerEdgeFToE2d = zeros(2, mesh2d.InnerEdge.Ne);
InnerEdgeFToE3d = zeros(2, mesh3d.InnerEdge.Ne);
LYCor = sort(mesh2d.BoundaryEdge.yb(:,LOrder2d));
RYCor = sort(mesh2d.BoundaryEdge.yb(:,ROrder2d));
for i = 1:size(LYCor,2)
    order = find(all(bsxfun(@eq,RYCor,LYCor(:,i)),1));
    ele1 = mesh2d.BoundaryEdge.FToE(1,LOrder2d(i));
    ele2 = mesh2d.BoundaryEdge.FToE(1,ROrder2d(order));
    InnerEdgeFToE2d(1,i) = ele1;
    InnerEdgeFToE2d(2,i) = ele2;
    for j = 1:mesh3d.Nz
        InnerEdgeFToE3d(1,(i-1)*mesh3d.Nz + j) = (ele1-1)*mesh3d.Nz + j;
        InnerEdgeFToE3d(2,(i-1)*mesh3d.Nz + j) = (ele2-1)*mesh3d.Nz + j;
    end
end
InnerEdgeFToE3d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.FToE;
InnerEdgeFToE2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToE;
mesh3d.InnerEdge.FToE = InnerEdgeFToE3d;
mesh2d.InnerEdge.FToE = InnerEdgeFToE2d;
mesh2d.BoundaryEdge.FToE(:,[LOrder2d, ROrder2d]) = [];
mesh3d.BoundaryEdge.FToE(:,[LOrder3d, ROrder3d]) = [];
%nx, ny and nz part
InnerEdgenx3d = zeros(size(mesh3d.InnerEdge.nx,1),mesh3d.InnerEdge.Ne);
InnerEdgeny3d = zeros(size(mesh3d.InnerEdge.ny,1),mesh3d.InnerEdge.Ne);
InnerEdgenz3d = zeros(size(mesh3d.InnerEdge.nz,1),mesh3d.InnerEdge.Ne);
InnerEdgenx3d(:,1:numel(LOrder3d)) = mesh3d.BoundaryEdge.nx(:,LOrder3d);
InnerEdgenx3d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.nx;
InnerEdgeny3d(:,1:numel(LOrder3d)) = mesh3d.BoundaryEdge.ny(:,LOrder3d);
InnerEdgeny3d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.ny;
InnerEdgenz3d(:,1:numel(LOrder3d)) = mesh3d.BoundaryEdge.nz(:,LOrder3d);
InnerEdgenz3d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.nz;

InnerEdgenx2d = zeros(size(mesh2d.InnerEdge.nx,1),mesh2d.InnerEdge.Ne);
InnerEdgeny2d = zeros(size(mesh2d.InnerEdge.ny,1),mesh2d.InnerEdge.Ne);
InnerEdgenx2d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.nx(:,LOrder2d);
InnerEdgenx2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.nx;
InnerEdgeny2d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.ny(:,LOrder2d);
InnerEdgeny2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.ny;

mesh2d.InnerEdge.nx = InnerEdgenx2d; mesh2d.InnerEdge.ny = InnerEdgeny2d;
mesh3d.InnerEdge.nx = InnerEdgenx3d; mesh3d.InnerEdge.ny = InnerEdgeny3d;
mesh3d.InnerEdge.nz = InnerEdgenz3d;
mesh2d.BoundaryEdge.nx(:,[LOrder2d, ROrder2d]) = [];
mesh2d.BoundaryEdge.ny(:,[LOrder2d, ROrder2d]) = [];
mesh3d.BoundaryEdge.nx(:,[LOrder3d, ROrder3d]) = [];
mesh3d.BoundaryEdge.ny(:,[LOrder3d, ROrder3d]) = [];
mesh3d.BoundaryEdge.nz(:,[LOrder3d, ROrder3d]) = [];
%FToN1 and FToN2 part
InnerEdgeFToN13d = zeros(size(mesh3d.InnerEdge.FToN1,1),mesh3d.InnerEdge.Ne);
InnerEdgeFToN23d = zeros(size(mesh3d.InnerEdge.FToN2,1),mesh3d.InnerEdge.Ne);
InnerEdgeFToN12d = zeros(size(mesh2d.InnerEdge.FToN1,1),mesh2d.InnerEdge.Ne);
InnerEdgeFToN22d = zeros(size(mesh2d.InnerEdge.FToN2,1),mesh2d.InnerEdge.Ne);
InnerEdgeFToN12d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.FToN1(:,LOrder2d);
InnerEdgeFToN22d(:,1:numel(ROrder2d)) = flip(mesh2d.BoundaryEdge.FToN2(:,ROrder2d));
InnerEdgeFToN13d(:,1:numel(LOrder3d)) = mesh3d.BoundaryEdge.FToN1(:,LOrder3d);
InnerEdgeFToN23d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.FToN2(:,ROrder3d);
InnerEdgeFToN12d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToN1;
InnerEdgeFToN22d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToN2;
InnerEdgeFToN13d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.FToN1;
InnerEdgeFToN23d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.FToN2;
mesh2d.InnerEdge.FToN1 = InnerEdgeFToN12d;
mesh2d.InnerEdge.FToN2 = InnerEdgeFToN22d;
mesh3d.InnerEdge.FToN1 = InnerEdgeFToN13d;
mesh3d.InnerEdge.FToN2 = InnerEdgeFToN23d;
mesh2d.BoundaryEdge.FToN1(:,[LOrder2d, ROrder2d]) = [];
mesh2d.BoundaryEdge.FToN2(:,[LOrder2d, ROrder2d]) = [];
mesh3d.BoundaryEdge.FToN1(:,[LOrder3d, ROrder3d]) = [];
mesh3d.BoundaryEdge.FToN2(:,[LOrder3d, ROrder3d]) = [];

%Js part
InnerEdgeJs3d = zeros(size(mesh3d.InnerEdge.Js,1),mesh3d.InnerEdge.Ne);
InnerEdgeJs2d = zeros(size(mesh2d.InnerEdge.Js,1),mesh2d.InnerEdge.Ne);
InnerEdgeJs2d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.Js(:,LOrder2d);
InnerEdgeJs3d(:,1:numel(LOrder3d)) = mesh3d.BoundaryEdge.Js(:,LOrder3d);
InnerEdgeJs2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.Js;
InnerEdgeJs3d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.Js;
mesh2d.InnerEdge.Js = InnerEdgeJs2d;
mesh3d.InnerEdge.Js = InnerEdgeJs3d;
mesh2d.BoundaryEdge.Js(:,[LOrder2d, ROrder2d]) = [];
mesh3d.BoundaryEdge.Js(:,[LOrder3d, ROrder3d]) = [];
end

