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
fprintf('\n----------------------------------------------------------\n')

fprintf('Compile GOTM source files.\n');

NdgCompileGOTMSourceCode();

fprintf('GOTM source files compiled.\n');

fprintf('\n----------------------------------------------------------\n')

fprintf('%s:: Compile mex files.\n', mfilename);
NdgConfigure(12)
fprintf('\n----------------------------------------------------------\n')

fprintf('%s:: Finish all the setup process.\n\n', mfilename);