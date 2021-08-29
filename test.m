% n = 100;
% As = cell(100,1);
% As{1} = sparse(rand(20,20));
% As{1} = As{1} + (As{1})';
% rhses = cell(100,1);
% rhses{1} = rand(20,1);
% for i = 2:100
%     As{i} = sparse(rand(20,20));
%     As{i} = As{i} + (As{i})';
%     rhses{i} = rand(20,1);
% end
% 
% pardiso_info = [];
% bSolverInitialized = 0;
% for i = 1:n
%     A = As{i};
%     rhs = rhses{i};
%     
%     A = A + eps*speye(size(A)); % diagonal must be full
%     A = tril(A);
%     % initialize
%     if ~bSolverInitialized
%         bSolverInitialized = 1;
%         pardiso_info = pardisoinit(-2,0);
%         % Analyze the matrix and compute a symbolic factorization.
%         pardiso_info = pardisoreorder(A, pardiso_info, false);
%     end
%     % Compute the numeric factorization
%     pardiso_info = pardisofactor(A, pardiso_info, false);
%     % Compute the solutions using the symbolic factorization
% %     [sol, obj.pardiso_info] = pardisosolve(A, rhs, pardiso_info, false);
%    [sol, ~] = pardisosolve(A, rhs, pardiso_info, false)
% end
% pardisofree(pardiso_info);
% clear pardiso_info;





Length = 1000;
n = 100;
As = cell(n,1);
As{1} = sparse(rand(Length,Length));
rhses = cell(n,1);
rhses{1} = rand(Length,1);
for i = 2:100
    As{i} = sparse(rand(Length,Length));
    rhses{i} = rand(Length,1);
end

pardiso_info = [];
bSolverInitialized = 0;
for i = 1:n
    A = As{i};
    rhs = rhses{i};
    
    A = A + eps*speye(size(A)); % diagonal must be full
    % initialize
    if ~bSolverInitialized
        bSolverInitialized = 1;
        pardiso_info = pardisoinit(11,0);
        % Analyze the matrix and compute a symbolic factorization.
        pardiso_info = pardisoreorder(A, pardiso_info, false);
    end
    % Compute the numeric factorization
    pardiso_info = pardisofactor(A, pardiso_info, false);
    % Compute the solutions using the symbolic factorization
%     [sol, obj.pardiso_info] = pardisosolve(A, rhs, pardiso_info, false);
   [sol, ~] = pardisosolve(A, rhs, pardiso_info, false)
end
pardisofree(pardiso_info);
clear pardiso_info;