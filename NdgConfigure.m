function NdgConfigure( varargin )
global COMPILER
% % initialize CFLAGS & LDFLAGS
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
    [path, 'JacobiP.c'],...
    [path, 'GradJacobiP.c']};
libfile = { [path, 'polylib.c'] };
FuncHandle(outPath, srcfile, libfile);
% NdgMesh
path = 'NdgMesh/@NdgMesh/private/';
srcfile = {[path, 'mxGetMeshIntegralValue.c'], ...
    [path, 'mxAssembleMeshConnection.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgMesh/@NdgMesh2d/private/';
srcfile = {[path, 'mxGetVertexWeightedData.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgMesh/@NdgExtendMesh3d/private/';
srcfile = {[path, 'mxGetMeshIntegralValue.c'],...
    [path,'mxExtend2dField.c'],...
    [path,'mxVerticalColumnIntegralField.c'],...
    [path, 'VerticalIntegralFromBottom.c'],...
    [path, 'mxGetVertexWeightedData.c']};
libfile = {'NdgMath/NdgMath.c',...
    'NdgMath/NdgSWE3D.c',...
    'NdgMath/NdgSWE.c'};
FuncHandle(path, srcfile, libfile);
% NdgEdge
path = 'NdgEdge/@NdgInnerEdge/private/';
srcfile = {[path, 'mxEvaluateStrongFromEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c'], ...
    [path, 'mxEvaluateStrongFormEdgeAlterRHS.c'], ...
    [path, 'mxEvaluateStrongFormEdgeCentralRHS.c'], ...
    [path, 'mxGetMeshIntegralValue.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgHaloEdge/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgSideEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c'], ...
    [path, 'mxGetMeshIntegralValue.c']};
% libfile = {};
libfile = {'NdgMath/NdgMath.c'};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgSideEdge3d/private/';
srcfile = {[path, 'mxVerticalColumnIntegralField.c'],...
    [path, 'mxVerticalRepmatFacialValue.c']};
libfile = {'NdgMath/NdgMath.c',...
    'NdgMath/NdgSWE.c',...
    'NdgMath/NdgSWE3D.c'};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgBottomInnerEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c'], ...
    [path, 'mxEvaluateStrongFormEdgeAlterRHS.c'], ...
    [path, 'mxEvaluateStrongFormEdgeCentralRHS.c'],...
    [path, 'mxGetMeshIntegralValue.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgBottomHaloEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c']};
libfile = {};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge/@NdgHaloEdge3d/private/';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c'], ...
    [path, 'mxEvaluateSurfValue.c']};
libfile = {'NdgMath/NdgMath.c'};
FuncHandle(path, srcfile, libfile);

path = 'NdgEdge\@NdgSurfaceHaloEdge3d\private\';
srcfile = {[path, 'mxEvaluateStrongFormEdgeRHS.c']};
libfile = {'NdgMath/NdgMath.c'};
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

path = 'NdgLimiter/NdgVertLimiter/@NdgVertLimiter3d/private/';
srcfile = {[path,'mxVertLimit3dNew.c'],...
    [path,'mxVertLimit3dWithLargeStencil.c']};
libfile = {'NdgMath/NdgMath.c',...
    'NdgMath/NdgSWE.c'};
FuncHandle(path, srcfile, libfile);

% SWE1d
path = 'Application/SWE/SWE1d/@SWEAbstract1d/private/';
AddIncludePath(path);
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateSurfaceValue1d.c'], ...
    [path, 'mxUpdateTimeInterval1d.c']};
FuncHandle(path, srcfile, libfile);
% Surface flux for SWE1d
path = 'Application/SWE/FaceFluxSolver/SWEFaceFluxSolver1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateSurfFlux1d.c']};
FuncHandle(path, srcfile, libfile);

% Surface flux for vertical momentum equation in the nonhydrostatic SWE1d
path = 'Application/SWE/FaceFluxSolver/SWENonhydroFaceFluxSolver1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateVerticalSurfFlux.c']};
FuncHandle(path, srcfile, libfile);

%Volume flux for SWE1d
path = 'Application/SWE/VolumeFluxSolver/SWEVolumeFluxSolver1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateFlux1d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/VolumeFluxSolver/SWEWDVolumeFluxSolver1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = { ...
    [path, 'mxEvaluateFlux1d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/VolumeFluxSolver/SWEPrebalanceVolumeFlux1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = {[path, 'mxEvaluateFlux1d.c']};
FuncHandle(path, srcfile, libfile);

%Volume flux of the vertical momentum equation for both the conventional and prebalanced SWE1d
path = 'Application/SWE/VolumeFluxSolver/SWENonhydroVolumeFluxSolver1d/private/';
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
srcfile = {[path, 'mxGetNonhydroVerticalVolumeFlux.c']};
FuncHandle(path, srcfile, libfile);

%Numerical flux for SWE1d
path = 'Application/SWE/NumFluxSolver/SWEHLLNumFluxSolver1d/private/';
srcfile = {[path, 'mxEvaluate.c']};
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWENonhydroHLLNumFluxSolver1d/private/';
srcfile = {[path, 'mxEvaluateUpwindNumFlux.c']};
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
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
AddIncludePath(path)
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c',...
    'NdgMath/NdgSWE.c',...
    'NdgMath/NdgMath.c'};
srcfile = { ...
    [path, 'mxImposeBoundaryCondition.c'], ...
    [path, 'mxHydrostaticReconstruction.c'], ...
    [path, 'mxUpdateTimeInterval2d.c'], ...
    [path, 'mxEvaluateBoundaryVertScope.c']};
FuncHandle(path, srcfile, libfile);

% Surface flux for SWE2d
path = 'Application/SWE/FaceFluxSolver/SWEFaceFluxSolver2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluateSurfFlux.c']};
FuncHandle(path, srcfile, libfile);

% Surface flux for vertical momentum equation in the nonhydrostatic SWE2d
path = 'Application/SWE/FaceFluxSolver/SWENonhydroFaceFluxSolver2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluateVerticalSurfFlux.c']};
FuncHandle(path, srcfile, libfile);

% Volume flux solver for the conventional SWE2d solver
path = 'Application/SWE/VolumeFluxSolver/SWEVolumeFluxSolver2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = {...
    [path, 'mxEvaluateFlux2d.c']};
FuncHandle(path, srcfile, libfile);

% Volume flux solver for the prebalanced SWE2d solver
path = 'Application/SWE/VolumeFluxSolver/SWEPrebalanceVolumeFlux2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = {...
    [path, 'mxEvaluateFlux2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/VolumeFluxSolver/SWEWDPrebalanceVolumeFluxSolver2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = {...
    [path, 'mxEvaluateFlux2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE2d/@SWEConventional2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluatePostFunc2d.c'], ...
    [path, 'mxEvaluateSourceTopography2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE2d/@SWEPreBlanaced2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = { ...
    [path, 'mxEvaluateSourceTopography2d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE2d/@SWEWDPreBlanaced2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = {[path, 'mxUpdateWDWetDryState.c'], ...
    [path, 'mxEvaluateSourceTopography2d.c']};
FuncHandle(path, srcfile, libfile);

% SWE numerical flux
path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver1d/private/';
AddIncludePath(path)
path = 'Application/SWE/NumFluxSolver/SWEHLLNumFluxSolver1d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);
% path = 'Application/SWE/NumFluxSolver/SWERoeNumFluxSolver1d/private/';
% srcfile = {[path, 'mxEvaluate.c']};
% FuncHandle(path, srcfile, libfile);

%

path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver2d/private/';
AddIncludePath(path);
path = 'Application/SWE/NumFluxSolver/SWEHLLNumFluxSolver2d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWEHLLCNumFluxSolver2d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);


path = 'Application/SWE/NumFluxSolver/SWELFNumFluxSolver2d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWEHLLCNumFluxSolver2d/private/';
srcfile = {[path, 'mxEvaluateNew.c']};
libfile = { 'NdgMath/NdgMath.c' ,...
    'NdgMath/NdgSWE.c'};
FuncHandle(path, srcfile, libfile);

% path = 'Application/SWE/NumFluxSolver/SWERoeNumFluxSolver2d/private/';
% srcfile = {[path, 'mxEvaluate.c']};
% FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver3d/private/';
AddIncludePath(path);
path = 'Application/SWE/NumFluxSolver/SWELFNumFluxSolver3d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWEAbstractNumFluxSolver3d/private/';
AddIncludePath(path);
path = 'Application/SWE/NumFluxSolver/SWEHLLNumFluxSolver3d/private/';
srcfile = {[path, 'mxEvaluate.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/NumFluxSolver/SWENonhydroHLLNumFluxSolver2d/private/';
srcfile = {[path, 'mxEvaluateUpwindNumFlux.c']};
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgNonhydrostaticSolver2d/private/';
srcfile = {[path,'mxAssembleWetDryInterface2d.c'],...
    [path,'mxEvaluateDownwindNumFlux.c'],...
    [path,'mxEvaluateSurfFlux.c'],...
    [path,'mxEvaluateUpwindNumFlux.c'],...
    [path,'mxGetElementFaceNormalDirectionVector.c'],...
    [path,'mxGetWetDryFace.c']};
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
% RemoveMexFile(path, srcfile, libfile);
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgAbstractNonhydrostaticSolver/private/';
srcfile = {[path,'mxAssembleElementBoundaryCondition.c'],...
    [path,'mxCalculatePenaltyParameter.c'],...
    [path,'mxGetPrimitiveVariableBoundaryEdgeFlux.c'],...
    [path,'mxGetPrimitiveVariableInnerEdgeFlux.c'],...
    [path,'mxSetBoundaryType.c']};
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
% RemoveMexFile(path, srcfile, libfile);
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver2d/private/';
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
libfile = {'Application/SWE/SWE1d/@SWEAbstract1d/private/mxSWE1d.c'};
% RemoveMexFile(path, srcfile, libfile);
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver1d/private/';
srcfile = {[path,'mxAssembleGlobalStiffMatrix.c'],...
    [path,'mxAssemblePointToCellInformation.c']};
libfile = { };
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/';
srcfile = {
    [path,'mxAssembleNonhydroRHS.c'],...
    [path,'mxUpdateConservativeFinalVelocity.c'],...
    [path,'mxCalculatePartialDerivativeUpdated.c'],...
    [path,'mxAssembleGlobalStiffMatrixNew.c'],...
    [path,'mxAssembleFinalGlobalStiffMatrix.c'],...
    [path,'mxFirstOrderDerivAboutNohydroPressInHorizontal.c'],...
    [path,'mxFirstOrderDerivAboutNohydroPressInVert.c'],...
    [path,'mxUpdateVerticalFinalVelocity.c'],...
    [path,'mxImposeBottomBoundaryCondition.c'],...
    [path,'mxImposeSideBoundaryCondition.c'],...
    [path,'mxCalculateBottomVerticalVelocity.c'],...
    [path,'mxAssemblePositiveDefiniteStiffMatrix.c'],...
    [path,'mxCalculateNonhydroRHS.c'],...
    [path,'mxGetInnerEdgeTopologyRelation.c'],...
    [path,'mxUpdateVerticalFinalVelocityIntegralForm.c']};
libfile = { [path, 'SWENonhydrostatic3d.c'],...
    'NdgMath/NdgMath.c',...
    'NdgMath/NdgSWE.c',...
    'NdgMath/NdgMemory.c',...
    'NdgMath/NdgSWE3D.c'};
FuncHandle(path, srcfile, libfile);

% path = 'Application\SWE\SWE3d\UnitTest\@EllipticProblemInVertical3d\private\';
% srcfile = {[path,'mxAssembleGlobalStiffMatrixWithSurfaceBCsImposed.c'],...
%     [path,'mxAssembleGlobalStiffMatrixWithBottomBCsImposed.c']};
% libfile = {'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/SWENonhydrostatic3d.c',...
%     'NdgMath/NdgMath.c',...
%     'NdgMath/NdgSWE.c',...
%     'NdgMath/NdgMemory.c',...
%     'NdgMath/NdgSWE3D.c'};
% FuncHandle(path, srcfile, libfile);
% 
% path = 'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/Unittest/NonhydrostaticSolver3dTest/@EllipticProblemMatrixAssembleTest3dNew/private/';
% srcfile = {[path,'mxAssembleGlobalStiffMatrixWithBCsImposed.c'],...
%     [path,'mxAssembleGlobalStiffMatrixWithSurfaceBCsImposed.c'],...
%     [path,'mxAssembleGlobalStiffMatrixWithBottomBCsImposed.c']};
% libfile = {'NdgPhys/NdgMatSolver/NdgNonhydrostaticSolver/@NdgQuadratureFreeNonhydrostaticSolver3d/private/SWENonhydrostatic3d.c',...
%     'NdgMath/NdgMath.c',...
%     'NdgMath/NdgSWE.c',...
%     'NdgMath/NdgMemory.c',...
%     'NdgMath/NdgSWE3D.c'};
% FuncHandle(path, srcfile, libfile);

% SWE2d

path = 'Application/SWE/VolumeFluxSolver/SWENonhydroVolumeFluxSolver2d/private/';
libfile = {'Application/SWE/SWE2d/@SWEAbstract2d/private/mxSWE2d.c'};
srcfile = {[path, 'mxGetNonhydroVerticalVolumeFlux.c']};
FuncHandle(path, srcfile, libfile);

% SWE3d
path = 'Application/SWE/SWE3d/@SWEAbstract3d/private/';
AddIncludePath(path);
libfile = {'Application/SWE/SWE3d/@SWEAbstract3d/private/mxSWE3d.c',...
    'NdgMath/NdgMath.c',...
    'NdgMath/NdgMemory.c'};
srcfile = { ...
    [path, 'mxUpdateTimeInterval3d.c'],...
    [path, 'mxUpdateNonhydroTimeInterval3d.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE3d/SWE3dVerticalVelocitySolver/private/';
AddIncludePath(path);
libfile = {'NdgMath/NdgMath.c',...
    'NdgMath/NdgSWE.c',...
    'NdgMath/NdgSWE3D.c',...
    'NdgMath/NdgMemory.c'};
srcfile = { ...
    [path, 'mxCalculateVerticalVelocity.c'],...
    [path, 'mxCalculateVerticalVelocityUpdate.c'],...
    [path, 'mxCalculateVerticalVelocityVersionOne.c'],...
    [path, 'mxCalculateVerticalVelocityUpdateVersionTwo.c'],...
    [path, 'mxCalculateVerticalVelocityUpdateVersionFour.c'],...
    [path, 'mxCalculateVerticalVelocityVersionSix.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE3d/@SWEQuadFreeStrongFormPCESolver2d/private/';
libfile = { 'NdgMath/NdgMath.c' ,...
    'NdgMath/NdgSWE.c',...
    'NdgMath/NdgSWE3D.c',...
    'NdgMath/NdgMemory.c'};
srcfile = { ...
    [path, 'mxEvaluatePCERHS.c'],...
    [path, 'mxEvaluatePCERHSUpdated.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE3d/@SWEAbstract3d/private/';
AddIncludePath(path);
libfile = {'Application/SWE/SWE3d/@SWEAbstract3d/private/mxSWE3d.c'};
srcfile = { ...
    [path, 'mxUpdateTimeInterval3d.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/@NdgPhysMat/private/';
AddIncludePath(path);
libfile = {'NdgMath/NdgMath.c'};
srcfile = { ...
    [path, 'mxAssembleSystemRHS.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgDiffSolver/@AbstractDiffSolver/private/';
AddIncludePath(path);
libfile = {};
srcfile = { ...
    [path, 'mxEvaluateSurfValue.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgDiffSolver/@NdgVertDiffSolver/private/';
AddIncludePath(path);
libfile = {'NdgMath/NdgMath.c',...
    'NdgMath/NdgMemory.c'};
srcfile = { ...
    [path, 'mxUpdateImplicitRHS.c']};
FuncHandle(path, srcfile, libfile);

% For different platforms, the path to mkl.h should be given explicitly
% path = 'NdgPhys/NdgMatSolver/NdgDiffSolver/@NdgVertDiffSolver/private/';
% mklpath = 'D:\Software\Intel\Install\mkl\2022.0.0\include\';
% AddIncludePath(mklpath);
% AddIncludePath(path);
% libfile = { 'NdgMath/NdgMemory.c',...
%     [path, 'mxImplicitVerticalEddyViscosity.c']};
% srcfile = { ...
%     [path, 'mxSparseVersionUpdateImplicitRHS.c']};
% FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgDiffSolver/@NdgSWEVertGOTMDiffSolver/private/';

switch computer('arch')
    case 'win64'
        Opath = [pwd,'/lib/GOTM/*.obj'];
    otherwise
        Opath = [pwd,'/lib/GOTM/*.o'];
end
file = dir(Opath);
libfile = [];
for i = 1:numel(file)
    libfile{i} = ...
        [pwd,'/lib/GOTM/',file(i).name];
end
libfile{numel(file)+1} = [path,'mxGOTM.c'];
libfile{numel(file)+2} = 'NdgMath/NdgMemory.c';
libfile{numel(file)+3} = 'NdgMath/NdgMath.c';
srcfile = {[path,'mxUpdateEddyViscosity.c'],...
    [path,'mxFilterData.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgDiffSolver/@NdgSWEHorizDiffSolver/private/';
libfile = { 'NdgMath/NdgSWE.c',...
    'NdgMath/NdgMath.c' ,...
    'NdgMath/NdgMemory.c' ,...
    [path,'HorizontalDiffusion.c']};
srcfile = { ...
    [path, 'mxEvaluateHorizontalDiffRHS.c']};
FuncHandle(path, srcfile, libfile);

path = 'NdgPhys/NdgMatSolver/NdgAdvSolver/NdgQuadFreeAdvSolver/private/';
AddIncludePath(path);
libfile = {'NdgMath/NdgMath.c',...
    'NdgMath/NdgSWE.c',...
    'NdgMath/NdgSWE3D.c',...
    'NdgMath/NdgMemory.c'};
srcfile = { ...
    [path, 'mxEvaluateQuadFreeStrongFormAdvSWE3dRHS.c']};
FuncHandle(path, srcfile, libfile);

path = 'Application/SWE/SWE3d/@SWEBaroclinic3d/private/';
libfile = {'NdgMath/NdgMath.c',...
    'NdgMath/NdgSWE3D.c',...
    'NdgMath/NdgSWE.c',...
    'NdgMath/NdgMemory.c'};
srcfile = { ...
    [path, 'mxCalculateBaroclinicTerm.c'],...
    [path, 'mxCalculateDensityField.c']};
FuncHandle(path, srcfile, libfile);

fprintf('/n%s:: Compiled all the mex files./n', mfilename);

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
global COMPFLAGS CFLAGS LDFLAGS COMPILER
switch computer('arch')
    case 'win64'
        for n = 1:numel(srcfile)
            if( isNeedCompile(outPath, srcfile{n}, libfile) )
                fprintf('\n%s:: Compiling source file - \n%s to %s.\n', ...
                    mfilename, srcfile{n}, outPath);
                fprintf('%s\nCOMPFLAGS=%s\nLDFLAGS=%s\n', COMPILER, COMPFLAGS, LDFLAGS);
                file = [srcfile(n), libfile{:}];
                mex('-g','-v','-largeArrayDims', COMPILER, COMPFLAGS, '-O', LDFLAGS, ...
                    file{:}, '-outdir', outPath);
%                 mex('-v','-largeArrayDims', COMPFLAGS, LDFLAGS, ...
%                     file{:}, '-outdir', outPath);
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
        %                 CFLAGS = [CFLAGS,' -fopenmp -DDG_THREADS=', ...
        %                     num2str(parallelThreadNum), ' '];
        %                 LDFLAGS = [LDFLAGS, ' -fopenmp '];
        %% if intel compiler adopted, /openmp command is used
        COMPFLAGS = [COMPFLAGS, ' /openmp -DDG_THREADS=', ...
            num2str(parallelThreadNum), ' '];
        LDFLAGS = [LDFLAGS, ' -LD:\Software\Intel\Install\mkl\2022.0.0\lib\intel64', ...
            ' /openmp -lmkl_intel_lp64_dll.lib -lmkl_core_dll.lib -lmkl_intel_thread_dll.lib'];
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
        COMPFLAGS = [COMPFLAGS, ' -I', path, ' '];
    case 'glnxa64'
        CFLAGS = [CFLAGS, ' -I', path, ' '];
end
end