clear, clc;
Nz = 1;
NLayer = [1 2 3 4 5 6];
for i = 1:numel(NLayer)
    Solver = StandingWaveInAClosedChannel(1,Nz,20,NLayer(i));
    Solver.matTimeSteppingLai;
    Solver.outputFile2d.mergeOutputResult;
    disp([Nz, NLayer(i)])
    Solver.NonhydroPostprocess;
end

Nz = 2;
NLayer = [1 2 3];
for i = 1:numel(NLayer)
    Solver = StandingWaveInAClosedChannel(1,Nz,20,NLayer(i));
    Solver.matTimeSteppingLai;
    Solver.outputFile2d.mergeOutputResult;
    disp([Nz, NLayer(i)])
    Solver.NonhydroPostprocess;
end

Nz = 3;
NLayer = [1 2];
for i = 1:numel(NLayer)
    Solver = StandingWaveInAClosedChannel(1,Nz,20,NLayer(i));
    Solver.matTimeSteppingLai;
    Solver.outputFile2d.mergeOutputResult;
    disp([Nz, NLayer(i)])
    Solver.NonhydroPostprocess;
end

