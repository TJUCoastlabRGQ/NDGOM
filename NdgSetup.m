function NdgSetup( varargin )
if (nargin == 1 && ( isa(varargin{1}, 'double')))
    DG_Threads = varargin{1};
elseif (nargin == 0)
    DG_Threads = 1;
else
    msgID = 'NdgConfigure:Unknown input.';
    msgtext = 'Unknown computer architecture.';
    throw( MException(msgID, msgtext) );
end
fprintf('\n----------------------------------------------------------\n')
fprintf('%s:: Setup environment path.\n', mfilename);
addpath( genpath( 'lib' ) );
addpath( genpath( 'NdgCell' ) );
addpath( genpath( 'NdgMesh') );
addpath( genpath( 'NdgEdge') );
addpath( genpath( 'NdgPhys') );
addpath( genpath( 'NdgNetCDF') );
addpath( genpath( 'NdgLimiter') );
addpath( genpath( 'Utilities') );
addpath( genpath( 'Application') );
addpath( genpath( 'NdgFilter') );
addpath( genpath( 'ConvergenceTest') );
addpath( genpath( 'Result'));
addpath( genpath( 'NdgMath'));
fprintf('\n----------------------------------------------------------\n')

fprintf('Compile GOTM source files.\n');

NdgCompileGOTMSourceCode();

fprintf('GOTM source files compiled.\n');

fprintf('\n----------------------------------------------------------\n')

fprintf('%s:: Compile mex files.\n', mfilename);
% If compiler is pointed, NdgConfigure is called as NdgConfigure(DG_Threads, Compiler);
NdgConfigure(DG_Threads);

fprintf('\n----------------------------------------------------------\n')

fprintf('%s:: Finish all the setup process.\n\n', mfilename);
end