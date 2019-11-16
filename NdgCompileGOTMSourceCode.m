function NdgCompileGOTMSourceCode()
if exist('.\lib\GOTM','dir')
    return;
else
    mkdir('.\lib\GOTM');
    
    copyfile('.\thirdParty\GOTM\src\turbulence\*.F90','.\lib\GOTM');
    copyfile('.\thirdParty\GOTM\src\util\*.F90','.\lib\GOTM');
    Files = {'.\lib\GOTM\test_time.F90','.\lib\GOTM\test_eqstate.F90',...
        '.\lib\GOTM\ode_solvers_template.F90'};
    DeleteFiles(Files);
    
    IncludePath = '-I.\thirdParty\GOTM\include';
    Outpath = '.\lib\GOTM';
    Files = {'.\lib\GOTM\turbulence.F90','.\lib\GOTM\util.F90',...
        '.\lib\GOTM\gotm_version.F90','.\lib\GOTM\eqstate.F90',...
        '.\lib\GOTM\tridiagonal.F90','.\lib\GOTM\time.F90'};
    
    CompileF90File(Files,IncludePath, Outpath);
    Files = dir('.\lib\GOTM\*.F90');
    FilesName = cell(1,numel(Files));
    for i = 1:numel(Files)
        FilesName{1,i} = strcat('.\lib\GOTM\',Files(i).name);
    end
    CompileF90File(FilesName,IncludePath, Outpath);
    DeleteFiles(FilesName);
end

end

function DeleteFiles(file)
for i = 1:numel(file)
    delete(file{i});
end
end

function CompileF90File(file, Ipath, Outpath)
for i = 1:numel(file)
    mex('-c','-v','-g','-largeArrayDims',file{i},Ipath,'-outdir', Outpath);
end
end