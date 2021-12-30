function matTimeExlicitStepping222(obj)
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys = obj.fphys;
%> allocate space for the rhs to be stored
ExplicitRHS1d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3);
DiffusionCoefficient = obj.miu * ones(size(obj.meshUnion(1).x));
dt = obj.dt;
% visual = Visual2d( obj.mesh2d );
visual = makeVisualizationFromNdgPhys( obj );
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys = fphys;
    tloc = time + c(1)*dt;
    obj.matUpdateExternalField( tloc );
    [ExplicitRHS1d(:,:,1)] = ...
        -1 * matCalculateExplicitRHSTerm(obj, fphys, tloc, DiffusionCoefficient);
    for intRK = 1:2
        tloc = time + c( intRK+1 ) * dt;
        obj.matUpdateExternalField( tloc );
        fphys{1}(:,:,1) = Tempfphys{1}(:,:,1) + dt * EXa(intRK+1,1) * ExplicitRHS1d(:,:,1) + dt * EXa(intRK+1,2) * ExplicitRHS1d(:,:,2);
        [ExplicitRHS1d(:,:,intRK + 1)] =  -1 * matCalculateExplicitRHSTerm(obj, fphys, tloc, DiffusionCoefficient);
    end
    fphys{1}(:,:,1) = Tempfphys{1}(:,:,1) + dt * EXb(1) * ExplicitRHS1d(:,:,1) + dt * EXb(2) * ExplicitRHS1d(:,:,2) + ...
        dt * EXb(3) * ExplicitRHS1d(:,:,3);
    
    
    ExplicitRHS1d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3);
    time = time + dt;
    display(time);
    %     obj.matUpdateOutputResult( time, fphys );
    timeRatio = time / ftime;
    visual.drawResult( fphys{1}(:, :, 1) );
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys = fphys;
% obj.matUpdateFinalResult( time, fphys );
% obj.outputFile.closeOutputFile();


% Nmesh = obj.Nmesh;
% [rk4a, rk4b, rk4c] = GetRKParamter();
% 
% time = obj.getOption('startTime');
% ftime = obj.getOption('finalTime');
% resQ = cell( obj.Nmesh, 1 );
% for n = 1:obj.Nmesh
%     resQ{n} = zeros( obj.meshUnion(n).cell.Np, obj.meshUnion(n).K, 1 );
% end
% fphys = obj.fphys;
% % init limiter and output file
% visual = makeVisualizationFromNdgPhys( obj );
% hwait = waitbar(0,'Runing MatSolver....');
% while( time < ftime )
%     if( time + dt > ftime )
%         dt = ftime - time;
%     end
%     
%     for intRK = 1:3
%         tloc = time + rk4c(intRK) * dt;
%         frhs = -1 * matCalculateExplicitRHSTerm(obj, fphys, tloc, DiffusionCoefficient);
%         
%         for n = 1:Nmesh
%             resQ{n} = rk4a(intRK)*resQ{n} + dt*frhs;
%             fphys{n}(:,:, 1) ...
%                 = fphys{n}(:,:, 1) + rk4b(intRK)*resQ{n};
%         end
%         
%     end
%     visual.drawResult( fphys{1}(:, :, 1) );   
%     time = time + dt;
%     timeRatio = time / ftime;
%     waitbar( timeRatio, hwait, ...
%         ['Runing MatSolver ', num2str( timeRatio ), '....']);
% end
% obj.fphys = fphys;

end

% function [rk4a, rk4b, rk4c] = GetRKParamter()
% rk4a = [ 0.0, 0.0, -1 ];
% rk4b = [ 0.5, 0.5, 1/3 ];
% rk4c = [ 0.0, 0.5, 1 ];
% end


function RHS1d = matCalculateExplicitRHSTerm(obj, fphys, tloc, DiffusionCoefficient)
meshUnion = obj.meshUnion(1);
%> Calculation of the penalty parameter
Tau = zeros(1, meshUnion.K);
Tau = matCalculatePenaltyParameter(meshUnion, DiffusionCoefficient, Tau);
qx = matCalculateAuxialaryVariable( obj, meshUnion, fphys );
RHS1d = matCalculateRHS1d(obj, meshUnion, Tau, qx, fphys, tloc);
end

function qx = matCalculateAuxialaryVariable( obj, meshUnion, fphys )
edge = meshUnion.InnerEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( fphys );
[ fluxM ] = edge.nx .* fm(:,:,1);
[ fluxP ] = edge.nx .* fp(:,:,1);
[ fluxS ] = 0.5 * edge.nx .* (fm(:,:,1) + fp(:,:,1));
[ qx ] = edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );

edge = meshUnion.BoundaryEdge;
[ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
[ fluxM ] = edge.nx .* fm(:,:,1);
fluxS = zeros(1,2);
fluxS(1) = edge.nx(1) * fm(1);%Newmann boundary
% fluxS(2) = edge.nx(2) * 1/sqrt(4*tloc+1)*exp(-(1-0.5).^2/miu/(4*tloc+1));
fluxS(2) = edge.nx(2) * obj.DirichExact(2);

[ qx ] = qx + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );

%Volume Intergral
qx = ...
     -qx +  meshUnion.rx.*( meshUnion.cell.Dr * fphys{1}(:,:,1) );

qx = obj.miu .* qx;
end

function RHS1d = matCalculateRHS1d(obj, meshUnion, Tau, qx, fphys, tloc)
sigma = cell(1);
sigma{1} = obj.miu .* (meshUnion.rx .*(meshUnion.cell.Dr * fphys{1}(:,:,1)));
Qx = cell(1);
Qx{1} = qx;


edge = meshUnion.InnerEdge;
[sigmafm, sigmafp] = edge.matEvaluateSurfValue(sigma);
[ fm, fp ] = edge.matEvaluateSurfValue( fphys );
[ Qxfm, Qxfp ] = edge.matEvaluateSurfValue( Qx );
[ fluxM ] = edge.nx .* Qxfm(:,:,1);
[ fluxP ] = edge.nx .* Qxfp(:,:,1);
[ fluxS ] = 0.5 * edge.nx .* (sigmafm(:,:,1) + sigmafp(:,:,1)) - edge.nx .* (fm(:,:,1) - fp(:,:,1)) * Tau(1);   %Here the directional vector should be included, where to include Tau
[ RHS1d ] = edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );


edge = meshUnion.BoundaryEdge;
[sigmafm, ~] = edge.matEvaluateSurfValue(sigma);
[ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
[ Qxfm, ~ ] = edge.matEvaluateSurfValue( Qx );
[ fluxM ] = edge.nx .* Qxfm(:,:,1);
% miu = 0.01; %this parameter should equal to the setted value exactly
% fluxS(1) = miu*1/sqrt(4*tloc+1)*(-2)*(0-0.5)/miu/(4*tloc+1)*exp(-(0-0.5).^2/miu/(4*tloc+1))*(-1);
fluxS(1) = obj.NewmannExact(1)*(-1);
% fluxS(2) = edge.nx(2) * sigmafm(2) - Tau(1) * (fm(2) - 1/sqrt(4*tloc+1)*exp(-(1-0.5).^2/miu/(4*tloc+1)));
fluxS(2) = edge.nx(2) * sigmafm(2) - Tau(1) * (fm(2) - obj.DirichExact(2));

[ RHS1d ] = RHS1d + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );


RHS1d = ...
     RHS1d -  meshUnion.rx.*( meshUnion.cell.Dr * qx );
end

function Tau = matCalculatePenaltyParameter(mesh, DiffusionCoefficient, Tau)
%> @brief Evaluating the penalty parameter used to penalize the jump between adjacet cell used in IPDG for second order operator
%>@detail In this version, the Interior Penalty Discontinuous Galerkin(IPDG) method is used to treat
%> the second order diffusion operator. To do so, the penalty parameter is calculated according to
%> [1] Shahbazi K. An explicit expression for the penalty parameter of the interior penalty method[J].
%> Journal of Computational Physics, 2005, 205(2): 401-407.
%> [2] Pestiaux A. Parameterization of subgrid-scale processes in finite element sea ice-ocean models[D].
%> UCL-Universit¨¦ Catholique de Louvain, 2015. pg:28.
%> The formula is '$\tau=\frac{(D_p+1)(D_p+d)}{d}\frac{n_0}{2}\frac{A}{V}\miu$'
%> @param[in] mesh3d The three-dimensional mesh object
%> @param[in] mesh2d The two-dimensional mesh object
%> @param[in] Tau The pre-allocated penalty parameter
%> @param[in] DiffusionCoefficient The scalar diffusion parameter
%> @param[out] Tau The calculated penalty parameter with size ( Nz + 1 ) * K2d
P = mesh.cell.N;
%> for prisms, number of faces is 5
n0 = 2;
%> here Nz stands for ratio between area of surface and volume of the studied cell
% Nz = mesh3d.Nz;
for i = 1:numel(Tau)
    Tau(i) = (P+1)*(P+1)/1*n0/2*1/mesh.LAV(1)*DiffusionCoefficient(1);
end
% for i = 1:mesh2d.K
%     %> The surface most face for each column
%     Tau(1,i) = 10*(P+1)*(P+1)/1*n0/2*Nz*max(DiffusionCoefficient(LeftEidM, (i-1)*Nz+1));
%     for j = 2:Nz
%         Tau(j,i) = 10*(P+1)*(P+3)/3*n0/2*Nz*max(max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1)),...
%             max(DiffusionCoefficient(UpEidM, (i-1)*Nz+j)));
%     end
%     %> The bottom most face for each column
%     Tau(Nz+1,i) = 10*(P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz));
% end
end



function [Explicita, Implicita, Explicitb, Implicitb, Parameterc] = GetRKParamter()
GAMA = (2-sqrt(2))/2;
delta = 1-1/(2*GAMA);
Parameterc = [0 GAMA 1];
Explicita = [0 0 0;
    GAMA 0 0;
    delta 1-delta 0];
Implicita = [GAMA 0;
    (1-GAMA) GAMA];
Explicitb = [delta 1-delta 0];
Implicitb = [1-GAMA GAMA];
end