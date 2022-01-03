classdef EllipticProblemTest1d < SWEConventional1d
    
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
        
        NewmannExact
        
        NewmannRightExact
        
        DirichletExact
        
        DirichletLeftExact
    end
 
    
 methods
     function obj = EllipticProblemTest1d(N, M)
           obj = obj@SWEConventional1d();
           obj.matGetFunction;
           [ mesh ] = makeUniformMesh( N, M );       
           obj.initPhysFromOptions( mesh );
           
           x = mesh.x;
           obj.ExactSolution = eval(obj.Cexact);
           
           x = mesh.x(1);
%            obj.NewmannExact = 100;
%            obj.NewmannExact = -1 * eval(obj.DiffCexact);
           obj.NewmannExact = -1*eval(obj.DiffCexact);
           obj.DirichletLeftExact = eval(obj.Cexact);       
           
           x = mesh.x(end);
           obj.DirichletExact = eval(obj.Cexact); 
           obj.NewmannRightExact = 1*eval(obj.DiffCexact);
     end
     
     function EllipticProblemSolve(obj)
         [ obj.StiffMatrix, LRHS, RRHS ] = obj.matAssembleGlobalStiffMatrix;
         obj.RHS = obj.matAssembleRightHandSide;
         % The direction vector is considered in the initialization stage
         obj.RHS(:,1) = obj.RHS(:,1) + LRHS;
         % Minus or plus is considered when assemble RRHS
         obj.RHS(:,end) = obj.RHS(:,end) + RRHS;
         obj.SimulatedSolution = obj.StiffMatrix\obj.RHS(:);
         
%          figure;
%          hold on;
%          plot(obj.meshUnion.x(1:5:end), obj.ExactSolution(1:5:end),'ro','LineWidth',1.5);
%          plot(obj.meshUnion.x(:), obj.SimulatedSolution(:),'k','LineWidth',1.5);
%          legend({'Exact','Simulation'});
%          legend('boxoff');
% %          title('$u=x(1-x)e^{2x}$','interpreter','Latex');
% %          title('$u=x(1-x)$','interpreter','Latex');
%          title('$u=sin(-\frac{\pi}{2}x)$','interpreter','Latex');
%          box on;
%          set(gca, 'Linewidth',1.5, 'Fontsize',12);
     end     
 end
    
  methods(Access = protected)
      
     [ Matrix, Lrhs, Rrhs ] = matAssembleGlobalStiffMatrix(obj);
      
      rhs = matAssembleRightHandSide(obj);
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            fphys{1} = zeros(obj.meshUnion.cell.Np, obj.meshUnion.K, 3);
            %             fphys{1}(:,:,1) = 1/obj.miu*exp(-(obj.mesh3d(1).z+0.5).^2);
        end
        
        function matGetFunction(obj)
            syms x;
            obj.Cexact = sin(-pi/2*x);
%             obj.Cexact = x*(1-x);
%             obj.Cexact = x*(1-x)*exp(2*x);
            obj.DiffCexact = diff(obj.Cexact, x);
            obj.SecondDiffCexact = diff(obj.DiffCexact, x);
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
xlim = [-1, 1];
bcType = [enumBoundaryCondition.Clamped, enumBoundaryCondition.SlipWall];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end