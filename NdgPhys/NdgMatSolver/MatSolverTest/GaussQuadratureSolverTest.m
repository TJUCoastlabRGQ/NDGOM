%> @brief Unit test for the Gauss Quadrature based solver classes.
%
%> @code
%>     addpath(pwd);
%>     results = runtests('NdgCell/test/StdCellTest.m');
%> @endcode
% ======================================================================
%> This class is part of the NDG-FEM software. 
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef GaussQuadratureSolverTest < matlab.unittest.TestCase
    
    properties(MethodSetupParameter)
        %> test cell types
%         type = {...
%             enumStdCell.Tri, ...
%             enumStdCell.Quad
%             }
        %> test cell orders
%         N = { 2 }
%         Nz = { 2 }
%         M = { 5 }
%         Mz = { 2 }
    end
    
    properties(Constant)
        %> tolerance
        tol = 1e-10;
    end
    
    properties
        %> cell object
        Solver
    end
    
    methods(TestMethodSetup)
        %> get the StdCell object
        function set_test_Solver(test) %(test, N, Nz, M, Mz)
            test.Solver = getSolver;
        end% func
    end
    
    methods(Test, ParameterCombination = 'sequential')
        function test_interpolation_matrix(test)
            test.verifyEqual(sum( test.Solver.advectionSolver.SBFVfq{1}, 2 ), ones( size(test.Solver.advectionSolver.SBFVfq{1},1) ,1 ), 'AbsTol', test.tol);
            test.verifyEqual(sum( test.Solver.advectionSolver.BOTFVfq{1}, 2 ), ones( size(test.Solver.advectionSolver.BOTFVfq{1},1) ,1 ), 'AbsTol', test.tol);
            test.verifyEqual(sum( test.Solver.advectionSolver.BBFVfq{1}, 2 ), ones( size(test.Solver.advectionSolver.BBFVfq{1},1) ,1 ), 'AbsTol', test.tol);
            test.verifyEqual(sum( test.Solver.advectionSolver.IEFVfq{1}, 2 ), ones( size(test.Solver.advectionSolver.IEFVfq{1},1) ,1 ), 'AbsTol', test.tol);
            test.verifyEqual(sum( test.Solver.advectionSolver.BEFVfq{1}, 2 ), ones( size(test.Solver.advectionSolver.BEFVfq{1},1) ,1 ), 'AbsTol', test.tol);
        end
        
    end% methods
end% classdef

function Solver = getSolver
     Solver = StandingWaveInAClosedChannel( 2, 2, 5, 2);
end