function NdgConfigure( varargin )
global COMPILER CFLAGS
% initialize CFLAGS & LDFLAGS
configureCompilerSetting();

if (nargin == 0)
    FuncHandle = @CompileMexFile;
elseif (nargin == 1) && ( isa(varargin{1}, 'double') )
    configureParallelSetting(varargin{1});
    FuncHandle = @CompileMexFile;
elseif (nargin == 1) && ( isa(varargin{1}, 'char') )
    if ( strcmp(varargin{1}, '-clear') )
        FuncHandle = @RemoveMexFile;
    end
elseif (nargin == 2)
    COMPILER = varargin{1};
    configureParallelSetting( varargin{2} )
    FuncHandle = @CompileMexFile;
else
    msgID = 'NdgConfigure:Unknown input.';
    msgtext = 'Unknown computer architecture.';
    throw( MException(msgID, msgtext) );
end

% set output path
outPath = 'lib/';
% Polylib
path = 'thirdParty/Polylib/';
srcfile = {[path, 'zwglj.c'], ...
    [path, 'zwgl.c'], ...
    [path, 'JacobiP.c'], ...
    [path, 'GradJacobiP.c'] };
libfile = { [path, 'polylib.c'] };
FuncHandle(outPath, srcfile, libfile);
% NdgMesh
path = 'NdgMesh/@NdgMesh/private/';
srcfile = {[path, 'mxGetMeshIntegralValue.c'], ...
    [path, 'mxAssembleMeshConnection.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgMesh\@NdgExtendMesh3d\private\';
srcfile = {[path, 'mxGetMeshIntegralValue.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);
% NdgEdge
path = 'NdgEdge/@NdgInnerEdge/private/';
srcfile = {[path, 'mxEvaluateStrongFromEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c'], ...
    [path, 'mxEvaluateStrongFormEdgeAlterRHS.c'], ...
    [path, 'mxEvaluateStrongFormEdgeCentralRHS.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgHaloEdge/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgSideEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFromEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgBottomInnerEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c'], ...
    [path, 'mxEvaluateStrongFormEdgeAlterRHS.c'], ...
    [path, 'mxEvaluateStrongFormEdgeCentralRHS.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgBottomHaloEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgHaloEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

% Limiter
path = 'NdgLimiter/NdgBJ/@NdgBJAbstract/private/';
srcfile = {[path, 'mxEvaluateVertAverage.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgLimiter/NdgBJ/@NdgBJ2d/private/';
srcfile = {[path, 'mxBJ2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgLimiter/NdgVertLimiter/@NdgVertLimiter/private/';
srcfile = {[path, 'mxEvaluateVertAverage.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgLimiter/NdgVertLimiter/@NdgVertLimiter2d/private/';
srcfile = {[path, 'mxVertLimit2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgLimiter/NdgVertLimiter/@NdgVertLimiter3d/private/';
srcfile = {[path,'mxVertLimit3d.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgFilter/@NdgAbstractFilter/private/';
srcfile = {[path,'mxFilterSolution.c']};
FuncHandle(path, srcfile, libfile);

% SWE1d
path = 'Application\SWE\SWE1d\@SWEAbstract1d\private\';
CFLAGS = [CFLAGS, ' -I', path, ' '];
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateSurfaceValue1d.c'], ...
    [path, 'mxUpdateTimeInterval1d.c']};
FuncHandle(path, srcfile, libfile);
% Surface flux for SWE1d
path = 'Application\SWE\FaceFluxSolver\SWEFaceFluxSolver1d\private\';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateSurfFlux1d.c']};
FuncHandle(path, srcfile, libfile);

% Surface flux for vertical momentum equation in the nonhydrostatic SWE1d
path = 'Application\SWE\FaceFluxSolver\SWENonhydroFaceFluxSolver1d\private\';
libfile = {'Application\SWE\SWE1d\@SWEAbstract1d\private\mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateVerticalSurfFlux.c']};
FuncHandle(path, srcfile, libfile);

%Volume flux for SWE1d
path = 'Application\SWE\VolumeFluxSolver\SWEVolumeFluxSolver1d\private\';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateFlux1d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application\SWE\VolumeFluxSolver\SWEWDVolumeFluxSolver1d\private\';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateFlux1d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application\SWE\VolumeFluxSolver\SWEPrebalanceVolumeFlux1d\private\';
libfile = {'Application\SWE\SWE1d\@SWEAbstract1d\private\mxSWE1d.c'};
srcfile = {[path, 'mxEvaluateFlux1d.c']};
FuncHandle(path, srcfile, libfile);

%Volume flux of the vertical momentum equation for both the conventional and prebalanced SWE1d
path = 'Application\SWE\VolumeFluxSolver\SWENonhydroVolumeFluxSolver1d\private\';
libfile = {'Application\SWE\SWE1d\@SWEAbstract1d\private\mxSWE1d.c'};
srcfile = {[path, 'mxGetNonhydroVerticalVolumeFlux.c']};
FuncHandle(path, srcfile, libfile);

%Numerical flux for SWE1d
path = 'Application\SWE\NumFluxSolver\SWEHLLNumFluxSolver1d\private\';
srcfile = {[path, 'mxEvaluate.c']};
libfile = {'Application\SWE\SWE1d\@SWEAbstract1d\private\mxSWE1d.c'};
FuncHandle(path, srcfile, libfile);

path = 'Application\SWE\NumFluxSolver\SWENonhydroHLLNumFluxSolver1d\private\';
srcfile = {[path, 'mxEvaluateUpwindNumFlux.c']};
libfile = {'Application\SWE\SWE1d\@SWEAbstract1d\private\mxSWE1d.c'};
FuncHandle(path, srcfile, libfile);

% Source function and post function for both conventional and prebalanced SWE1d
path = 'Application/SWE/SWE1d/@SWEConventional1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = {...
    [path, 'mxEvaluateSourceTopography1d.c'], ...
    [path, 'mxEvaluatePostFunc1d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE1d/@SWEPreBlanaced1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = {[path, 'mxEvaluateSourceTopography1d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE1d/@SWEWD1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = {[path, 'mxEvaluateSourceTopography1d.c']};
FuncHandle(path, srcfile, libfile);

% SWE2d
path = 'Application/SWE/SWE2d/@SWEAbstract2d/private/';
CFLAGS = [CFLAGS, ' -I', path, ' '];
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = { ...
    [path, 'mxImposeBoundaryCondition.c'], ...
    [path, 'mxHydrostaticReconstruction.c'], ...
    [path, 'mxUpdateTimeInterval2d.c']};
FuncHandle(path, srcfile, libfile);

% Surface flux for SWE2d
path = 'Application\SWE\FaceFluxSolver\SWEFaceFluxSolver2d\private\';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluateSurfFlux.c']};
FuncHandle(path, srcfile, libfile);

% Surface flux for vertical momentum equation in the nonhydrostatic SWE2d
path = 'Application\SWE\FaceFluxSolver\SWENonhydroFaceFluxSolver2d\private\';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluateVerticalSurfFlux.c']};
FuncHandle(path, srcfile, libfile);

% Volume flux solver for the conventional SWE2d solver
path = 'Application\SWE\VolumeFluxSolver\SWEVolumeFluxSolver2d\private\';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = {...
    [path, 'mxEvaluateFlux2d.c']};
FuncHandle(path, srcfile, libfile);

% Volume flux solver for the prebalanced SWE2d solver
path = 'Application\SWE\VolumeFluxSolver\SWEPrebalanceVolumeFlux2d\private\';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = {...
    [path, 'mxEvaluateFlux2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application\SWE\VolumeFluxSolver\SWEWDPrebalanceVolumeFluxSolver2d\private\';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = {...
    [path, 'mxEvaluateFlux2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE2d/@SWEConventional2d/private/';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluatePostFunc2d.c'], ...
    [path, 'mxEvaluateSourceTopography2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE2d/@SWEPreBlanaced2d/private/';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluateSourceTopography2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE2d/@SWEWDPreBlanaced2d/private/';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = {[path, 'mxUpdateWDWetDryState.c'], ...
    [path, 'mxEvaluateSourceTopography2d.c']};
FuncHandle(path, srcfile, libfile);

% SWE numerical flux
path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver1d/private/';
CFLAGS = [CFLAGS, ' -I', path, ' '];
path = 'Application/SWE/NumFluxSolver/SWEHLLNumFluxSolver1d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);
% path = 'Application/SWE/NumFluxSolver/SWERoeNumFluxSolver1d/private/';
% srcfile = {[path, 'mxEvaluate.c']};
% FuncHandle(path, srcfile, libfile);

%

path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver2d/private/';
CFLAGS = [CFLAGS, ' -I', path, ' '];
path = 'Application/SWE/NumFluxSolver/SWEHLLNumFluxSolver2d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);
path = 'Application/SWE/NumFluxSolver/SWELFNumFluxSolver2d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);
% path = 'Application/SWE/NumFluxSolver/SWERoeNumFluxSolver2d/private/';
% srcfile = {[path, 'mxEvaluate.c']};
% FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver3d/private/';
CFLAGS = [CFLAGS, ' -I', path, ' '];
path = 'Application/SWE/NumFluxSolver/SWELFNumFluxSolver3d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver3d/private/';
CFLAGS = [CFLAGS, ' -I', path, ' '];
path = 'Application/SWE/NumFluxSolver/SWEHLLNumFluxSolver3d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application\SWE\NumFluxSolver\SWENonhydroHLLNumFluxSolver2d\private\';
srcfile = {[path, 'mxEvaluateUpwindNumFlux.c']};
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgNonhydrostaticSolver2d/private/';
srcfile = {[path,'mxAssembleWetDryInterface2d.c'],...
    [path,'mxEvaluateDownwindNumFlux.c'],...
    [path,'mxEvaluateSurfFlux.c'],...
    [path,'mxEvaluateUpwindNumFlux.c'],...
    [path,'mxGetElementFaceNormalDirectionVector.c'],...
    [path,'mxGetWetDryFace.c']};
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
% RemoveMexFile(path, srcfile, libfile);
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgAbstractNonhydrostaticSolver/private/';
srcfile = {[path,'mxAssembleElementBoundaryCondition.c'],...
    [path,'mxCalculatePenaltyParameter.c'],...
    [path,'mxGetPrimitiveVariableBoundaryEdgeFlux.c'],...
    [path,'mxGetPrimitiveVariableInnerEdgeFlux.c'],...
    [path,'mxSetBoundaryType.c']};
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
% RemoveMexFile(path, srcfile, libfile);
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys\NdgMatSolver\NdgNonhydrostaticSolver\@NdgQuadratureFreeNonhydrostaticSolver2d\private\';
srcfile = {[path,'mxAssembleGlobalStiffMatrix.c'],...
    [path,'mxAssemblePointToCellInformation.c']};
libfile = { };
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgNonhydrostaticSolver1d/private/';
srcfile = {[path,'mxAssembleWetDryInterface1d.c'],...
    [path,'mxEvaluateDownwindNumFlux.c'],...
    [path,'mxEvaluateSurfFlux.c'],...
    [path,'mxEvaluateUpwindNumFlux.c'],...
    [path,'mxGetElementFaceNormalDirectionVector.c'],...
    [path,'mxGetWetDryFace.c']};
libfile = {'Application\SWE\SWE1d\@SWEAbstract1d\private\mxSWE1d.c'};
% RemoveMexFile(path, srcfile, libfile);
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys\NdgMatSolver\NdgNonhydrostaticSolver\@NdgQuadratureFreeNonhydrostaticSolver1d\private\';
srcfile = {[path,'mxAssembleGlobalStiffMatrix.c'],...
    [path,'mxAssemblePointToCellInformation.c']};
libfile = { };
FuncHandle(path, srcfile, libfile);

% SWE2d
path = 'Application/SWE/SWE3d/@SWEAbstract3d/private/';
CFLAGS = [CFLAGS, ' -I', path, ' '];
libfile = {'Application/SWE/SWE3d/@SWEAbstract3d/private/mxSWE3d.c'};
srcfile = { ...
    [path, 'mxUpdateTimeInterval3d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application\SWE\VolumeFluxSolver\SWENonhydroVolumeFluxSolver2d\private/';
libfile = {'Application\SWE\SWE2d\@SWEAbstract2d\private\mxSWE2d.c'};
srcfile = {[path, 'mxGetNonhydroVerticalVolumeFlux.c']};
FuncHandle(path, srcfile, libfile);

% SWE3d
path = 'Application/SWE/SWE3d/@SWEAbstract3d/private/';
CFLAGS = [CFLAGS, ' -I', path, ' '];
libfile = {'Application/SWE/SWE3d/@SWEAbstract3d/private/mxSWE3d.c'};
srcfile = { ...
    [path, 'mxUpdateTimeInterval3d.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys\NdgMatSolver\NdgDiffSolver\@AbstractDiffSolver\private\';
CFLAGS = [CFLAGS, ' -I', path, ' '];
libfile = {};
srcfile = { ...
    [path, 'mxEvaluateSurfValue.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys\NdgMatSolver\NdgDiffSolver\@NdgVertDiffSolver\private\';
CFLAGS = [CFLAGS, ' -I', path, ' '];
libfile = {};
srcfile = { ...
    [path, 'mxUpdateImplicitRHS.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys\NdgMatSolver\NdgDiffSolver\@NdgSWEVertGOTMDiffSolver\private\';
libfile = {['.\lib\GOTM\','*.obj'],...
    [path,'mxGOTM.c']};
srcfile = {[path,'mxUpdateEddyViscosity.c']};
FuncHandle(path, srcfile, libfile);
fprintf('\n%s:: Compiled all the mex files.\n', mfilename);

end

function RemoveMexFile(outPath, srcfile, libfile)
for n = 1:numel(srcfile)
    mexFile = getMexFileName(outPath, srcfile{n});
    file = dir(mexFile);
    if ~isempty(file)
        fprintf('\n%s:: Removing mex file %s in %s.\n', ...
            mfilename, mexFile, outPath);
        delete( fullfile(pwd, mexFile) );
    end
end
end

function CompileMexFile(outPath, srcfile, libfile)
global CFLAGS LDFLAGS COMPILER
for n = 1:numel(srcfile)
    if( isNeedCompile(outPath, srcfile{n}, libfile) )
        fprintf('\n%s:: Compiling source file - \n%s to %s.\n', ...
            mfilename, srcfile{n}, outPath);
        fprintf('%s\nCFLAGS=%s\nLDFLAGS=%s\n', COMPILER, CFLAGS, LDFLAGS);
        file = [srcfile(n), libfile{:}];
        mex('-v',  '-largeArrayDims', COMPILER, CFLAGS, '-O', LDFLAGS, ...
            file{:}, '-outdir', outPath);        
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
global CFLAGS LDFLAGS
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
        CFLAGS = [CFLAGS,' /openmp -DDG_THREADS=', ...
            num2str(parallelThreadNum), ' '];
        LDFLAGS = [LDFLAGS, ' /openmp '];
    case 'glnxa64'
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