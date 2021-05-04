function matCompileFile( obj )
global COMPILER
% initialize CFLAGS & LDFLAGS
configureCompilerSetting();

FuncHandle = @CompileMexFile;

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/Unittest/@NonhydroStaticSolver3dTest/private/';
srcfile = {[path,'mxCalculateBoundaryEdgeNumFluxTerm.c'],...
    [path,'mxCalculateInnerEdgeNumFluxTerm.c']};
libfile = {'NdgMath/NdgSWE.c',...
    'NdgMath/NdgSWE3D.c',...
    'NdgMath/NdgMath.c'};
FuncHandle(path, srcfile, libfile);

fprintf('\n%s:: Compiled all the mex files.\n', mfilename);

end

function CompileMexFile(outPath, srcfile, libfile)
global COMPFLAGS CFLAGS LDFLAGS COMPILER
switch computer('arch')
    case 'win64'
        for n = 1:numel(srcfile)
            if( isNeedCompile(outPath, srcfile{n}, libfile) )
                fprintf('\n%s:: Compiling source file - \n%s to %s.\n', ...
                    mfilename, srcfile{n}, outPath);
                fprintf('%s\nCOMPFLAGS=%s\nLDFLAGS=%s\n', COMPILER, COMPFLAGS, LDFLAGS);
                file = [srcfile(n), libfile{:}];
                mex('-g','-largeArrayDims', COMPILER, COMPFLAGS, '-O', LDFLAGS, ...
                    file{:}, '-outdir', outPath);
            end
        end
    otherwise
        for n = 1:numel(srcfile)
            if( isNeedCompile(outPath, srcfile{n}, libfile) )
                fprintf('\n%s:: Compiling source file - \n%s to %s.\n', ...
                    mfilename, srcfile{n}, outPath);
                fprintf('%s\nCFLAGS=%s\nLDFLAGS=%s\n', COMPILER, CFLAGS, LDFLAGS);
                file = [srcfile(n), libfile{:}];
                mex('-v','-largeArrayDims', COMPILER, CFLAGS, '-O', LDFLAGS, ...
                    file{:}, '-outdir', outPath);
            end
        end
end
end

function flag = isNeedCompile(outPath, srcFile, libSrdFile)
mexFile = getMexFileName(outPath, srcFile);
file = dir(mexFile);
if isempty(file)
    flag = true;
    return
end
mexDateNum = file.datenum;
file = dir(srcFile);
srcDateNum = file.datenum;
for n = 1:numel(libSrdFile)
    file = dir(libSrdFile{n});
    % find the lasted srcDateNum
    for m = 1:numel(file)
        srcDateNum = max(file(m).datenum, srcDateNum);
    end
end
flag = (srcDateNum > mexDateNum);
end

%> return mex file name based on the source file
function mexFile = getMexFileName(path, srcFile)
[ ~,name,~ ] = fileparts(srcFile);
switch computer('arch')
    case 'maci64'
        ext = '.mexmaci64';
    case 'win64'
        ext = '.mexw64';
    case 'glnxa64'
        ext = '.mexa64';
end
mexFile = [path, name, ext];
end% func

function configureParallelSetting(parallelThreadNum)
global CFLAGS  COMPFLAGS LDFLAGS
switch computer('arch')
    case 'maci64'
        CFLAGS = [CFLAGS, '-qopenmp -DDG_THREADS=', ...
            num2str(parallelThreadNum), ' '];
        LDFLAGS = [LDFLAGS, ' -lmwblas -liomp5 '];
    case 'win64'
        %% if gcc compiler adopted, -fopenmp command is used
        %         CFLAGS = [CFLAGS,' -fopenmp -DDG_THREADS=', ...
        %             num2str(parallelThreadNum), ' '];
        %         LDFLAGS = [LDFLAGS, ' -fopenmp '];
        %% if intel compiler adopted, /openmp command is used
        COMPFLAGS = [COMPFLAGS,' /openmp -DDG_THREADS=', ...
            num2str(parallelThreadNum), ' '];
        LDFLAGS = [LDFLAGS, ' /openmp '];
    case 'glnxa64'
        CFLAGS = [CFLAGS,' -fopenmp -DDG_THREADS=', ...
            num2str(parallelThreadNum)];
        LDFLAGS = [LDFLAGS, ' -fopenmp '];
end
end

function configureCompilerSetting()
global CFLAGS LDFLAGS
global COMPILER
CFLAGS = 'CFLAGS=$CFLAGS -std=c99 -Wall -DPROFILE -largeArrayDims';
LDFLAGS = 'LDFLAGS=$LDFLAGS';
COMPILER = [''];
if ( strcmp(computer('arch'), 'maci64') )
    COMPILER = 'CC=/opt/intel/composer_xe_2015.5.222/bin/intel64/icc';
end
end

function AddIncludePath(path)
global CFLAGS  COMPFLAGS
switch computer('arch')
    case 'maci64'
        CFLAGS = [CFLAGS, ' -I', path, ' '];
    case 'win64'
        COMPFLAGS = [COMPFLAGS, ' -I', path,  ' '];
    case 'glnxa64'
        CFLAGS = [CFLAGS, ' -I', path, ' '];
end
end