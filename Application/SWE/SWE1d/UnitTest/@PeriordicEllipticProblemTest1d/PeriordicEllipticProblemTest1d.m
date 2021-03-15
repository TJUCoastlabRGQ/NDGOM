classdef PeriordicEllipticProblemTest1d < SWEConventional1d
    
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 0.01
        %> gravity acceleration
        gra = 9.8
    end
    
    properties
        
        StiffMatrix
        
        RHS
        
        ExactSolution
        
        SimulatedSolution
        
        ExactRHS
        
    end
    
    properties
        
        Cexact
        
        DiffCexact
        
        SecondDiffCexact
        
    end
 
    
 methods
     function obj = PeriordicEllipticProblemTest1d(N, M)
           obj = obj@SWEConventional1d();
           obj.matGetFunction;
           [ mesh ] = makeUniformMesh( N, M );       
           obj.initPhysFromOptions( mesh );        
           x = mesh.x;
           obj.ExactSolution = eval(obj.Cexact);    
     end
     
     function EllipticProblemSolve(obj)
         [ obj.StiffMatrix ] = obj.matAssembleGlobalStiffMatrix;
         obj.RHS = obj.matAssembleRightHandSide;
         % The direction vector is considered in the initialization stage
%          obj.RHS(:,1) = obj.RHS(:,1) - LRHS;
         % Minus or plus is considered when assemble RRHS
%          obj.RHS(:,end) = obj.RHS(:,end) + RRHS;
         obj.SimulatedSolution = obj.StiffMatrix\obj.RHS(:);
         
         obj.ExactRHS = obj.StiffMatrix * obj.ExactSolution(:);
         figure;
         hold on;
         plot(obj.meshUnion.x(1:20:end), obj.ExactSolution(1:20:end),'ro');
         plot(obj.meshUnion.x(:), obj.SimulatedSolution(:),'k','LineWidth',1.5);
         legend({'Exact','Simulation'});
         legend('boxoff');
%          title('$u=x(1-x)e^{2x}$','interpreter','Latex');
%          title('$u=x(1-x)$','interpreter','Latex');
         title('$u=sin(-\frac{\pi}{2}x)$','interpreter','Latex');
         box on;
         set(gca, 'Linewidth',1.5, 'Fontsize',12);
     end     
 end
    
  methods(Access = protected)
      
     [ Matrix ] = matAssembleGlobalStiffMatrix(obj);
      
      rhs = matAssembleRightHandSide(obj);
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            fphys{1} = zeros(obj.meshUnion.cell.Np, obj.meshUnion.K, 3);
            %             fphys{1}(:,:,1) = 1/obj.miu*exp(-(obj.mesh3d(1).z+0.5).^2);
        end
        
        function matGetFunction(obj)
            syms x;
            obj.Cexact = sin(2*pi*x/20 + pi/2);
%             obj.Cexact = x*(1-x);
%             obj.Cexact = x*(1-x)*exp(2*x);
            obj.SecondDiffCexact = diff(diff(obj.Cexact, x),x);
        end        
        
        function [ option ] = setOption( obj, option )
            ftime = 40;
            outputIntervalNum = 5000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('EddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('ConstantEddyViscosityValue') = 0.01;
            option('limiterType') = enumLimiter.None;
            option('outputType') = enumOutputFile.NetCDF;
        end
    end    
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 10];
bcType = [enumBoundaryCondition.SlipWall, enumBoundaryCondition.SlipWall];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
% mesh.InnerEdge.Ne = mesh.InnerEdge.Ne + 1;
% mesh.BoundaryEdge.ftype = zeros(2,0);
% mesh.BoundaryEdge.Ne = 0;
% FToV = zeros(1, numel(mesh.InnerEdge.FToV) + 1 );
% FToV(1) = mesh.BoundaryEdge.FToV(1);
% % FToV(end) = mesh.BoundaryEdge.FToV(2);
% FToV(2:end) = mesh.InnerEdge.FToV;
% mesh.InnerEdge.FToV = FToV;
% mesh.BoundaryEdge.FToV = zeros(0,2);
% 
% FToE = zeros(2, mesh.InnerEdge.Ne);
% FToE(1) = mesh.K; FToE(2) = 1;
% LAV = zeros(1, mesh.InnerEdge.Ne);
% LAV(1) = mesh.BoundaryEdge.LAV(1);
% LAV(2:end) = mesh.InnerEdge.LAV;
% FToE(3:end) = mesh.InnerEdge.FToE(:);
% mesh.InnerEdge.FToE = FToE;
% mesh.InnerEdge.LAV = LAV;
% mesh.BoundaryEdge.FToE = zeros(2,0);
% mesh.BoundaryEdge.LAV = zeros(2,0);
% 
% FToF = zeros(2, mesh.InnerEdge.Ne);
% FToF(1) = 2; FToF(2) = 1;
% FToF(3:end) = mesh.InnerEdge.FToF(:);
% mesh.InnerEdge.FToF = FToF;
% mesh.BoundaryEdge.FToF = zeros(2,0);
% 
% FToN1 = zeros(1,mesh.InnerEdge.Ne);
% FToN2 = zeros(1,mesh.InnerEdge.Ne);
% FToN1(1) = mesh.BoundaryEdge.FToN1(2);
% FToN1(2:end) = mesh.InnerEdge.FToN1;
% mesh.InnerEdge.FToN1 = FToN1;
% FToN2(1) = mesh.BoundaryEdge.FToN1(1);
% FToN2(2:end) = mesh.InnerEdge.FToN2;
% mesh.InnerEdge.FToN2 = FToN2;
% mesh.BoundaryEdge.FToN1 = zeros(2,0);
% mesh.BoundaryEdge.FToN2 = zeros(2,0);
% 
% mesh.InnerEdge.nx = ones(1, mesh.InnerEdge.Ne);
% mesh.BoundaryEdge.nx = zeros(2,0);
% mesh.InnerEdge.ny = zeros(1, mesh.InnerEdge.Ne);
% mesh.BoundaryEdge.ny = zeros(0,2);
% mesh.InnerEdge.nz = zeros(1, mesh.InnerEdge.Ne);
% mesh.BoundaryEdge.nz = zeros(0,2);
% mesh.InnerEdge.Js = ones(1, mesh.InnerEdge.Ne);
% mesh.BoundaryEdge.Js = zeros(0,2);
% 
% mesh.EToE(1) = mesh.K;
% mesh.EToE(end) = 1;
% 
% mesh.EToF(1) = 2;
% mesh.EToF(end) = 1;
end