function MergeNcFile
ncfile = dir('*.nc');
Str = strsplit(ncfile(1).name,'.');
InputFile = strcat(Str{1},'.',Str{2},'.1.nc');
ncid1 = netcdf.open( InputFile, 'WRITE');
outputStep = numel(netcdf.getVar(ncid1,0));
for i = 2:numel(ncfile)
    InputFile = strcat(Str{1},'.',Str{2},'.',num2str(i),'.nc');
    ncid = netcdf.open( InputFile, 'WRITE');
    Time = netcdf.getVar(ncid,0);
    Info = ncinfo(InputFile);
    for n = 1:numel(Time)
        startInd = outputStep;
        countInd = 1;
        netcdf.putVar(ncid1, 0, startInd, countInd, Time(n));
        for m = 1:numel(Info.Variables) - 1
            field = netcdf.getVar(ncid,m);
            startInd = [ 0, 0, 0, outputStep ];
            countInd = [ Info.Variables( m+1 ).Size(1),  Info.Variables( m+1 ).Size(2),  Info.Variables( m+1 ).Size(3), 1 ];
            netcdf.putVar(ncid1, 1, startInd, countInd, field(:,:,:,n));            
        end
        outputStep = outputStep + 1;
    end
    netcdf.close(ncid);
    delete(InputFile);
end
netcdf.close(ncid1);
end