function runStandingWaveInAClosedBasin
close all;
clear, clc;
Nz = 1;
NL = [1 2 3 4 5 6];

[ ErrOrder1, TimeOrder1 ] = runSimulations( Nz, NL);
disp('=======Errors for order 1================')
disp(ErrOrder1)
disp('=======End errors for order 1================')
disp('=======Times for order 1================')
disp(TimeOrder1)
disp('=======End times for order 1================')

Nz = 2;
NL = [1 2 3];

[ ErrOrder2, TimeOrder2 ] = runSimulations( Nz, NL);
disp('=======Errors for order 2================')
disp(ErrOrder2)
disp('=======End errors for order 2================')
disp('=======Times for order 2================')
disp(TimeOrder2)
disp('=======End times for order 2================')

Nz = 3;
NL = [1 2];

[ ErrOrder3, TimeOrder3 ] = runSimulations( Nz, NL);
disp('=======Errors for order 3================')
disp(ErrOrder3)
disp('=======End errors for order 3================')
disp('=======Times for order 3================')
disp(TimeOrder3)
disp('=======End times for order 3================')

end

function [Err, Time] = runSimulations( Nz, NL )
Time = zeros(numel(NL),1);
Err = zeros(numel(NL),1);
for i = 1:numel(NL)
    Solver = StandingWaveInAClosedChannel(1, Nz, 40, NL(i));
    tic;
    Solver.matTimeSteppingLai;
    Time(i) = toc;
    Solver.outputFile2d.mergeOutputResult;
    PostProcess = NdgPostProcess(Solver.meshUnion(1).mesh2d,strcat('Result/','StandingWaveInAClosedChannel','/2d/','StandingWaveInAClosedChannel'));
    Ntime = PostProcess.Nt;
    outputTime = ncread( PostProcess.outputFile{1}, 'time' );
    Eta = zeros( Ntime,1 );
    exactEta = zeros( Ntime,1 );
    x0 = 18;
    h = Solver.H0;
    a = Solver.A;
    c = sqrt( Solver.gra*Solver.Lambda/2/pi*tanh(2*pi*h/Solver.Lambda) );
    T = Solver.Lambda / c;
    for t = 1:Ntime
        %                 exactEta(t) = -obj.A * cos(2*pi*x0/obj.Lambda + 2*pi*outputTime(t) /obj.T);
        exactEta(t) = Solver.A * cos(2*pi/Solver.Lambda*x0)*cos(2*pi/T*outputTime(t));
        tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  x0, 0.2, x0, t )-Solver.H0;
        Eta(t) = tempdata(1);
        Err(i) = Err(i) + (exactEta(t)  - Eta(t))^2;
    end
    Err(i) = sqrt(Err(i)/Ntime);
    
end
end