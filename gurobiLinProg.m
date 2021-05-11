function [x, fval, exitflag] = linprog(f, A, b, Aeq, beq, lb, ub)
% A linear programming formalization using the Gurobi MATLAB interface
% This code was downloaded from Gruobi's website.

if nargin < 3
    error('linprog(f, A, b)')
end

if nargin > 7
    error('linprog(f, A, b, Aeq, beq, lb, ub)');
end

if ~isempty(A)
    n = size(A, 2);
elseif nargin > 4 && ~isempty(Aeq)
    n = size(Aeq, 2);
else
    error('No linear constraints specified')
end

if ~issparse(A)
    A = sparse(A);
end

if nargin > 3 && ~issparse(Aeq)
    Aeq = sparse(Aeq);
end


model.obj = f;

if nargin < 4
    model.A = A;
    model.rhs = b;
    model.sense = '<';
else
    model.A = [A; Aeq];
    model.rhs = [b; beq];
    model.sense = [repmat('<', size(A,1), 1); repmat('=', size(Aeq,1), 1)];
end

if nargin < 6
    model.lb = -inf(n,1);
else
    model.lb = lb;
end

if nargin == 7
   model.ub = ub;
end

params.outputflag = 0;
result = gurobi(model, params);


if strcmp(result.status, 'OPTIMAL')
    exitflag = 1;
elseif strcmp(result.status, 'ITERATION_LIMIT')
    exitflag = 0;
elseif strcmp(result.status, 'INF_OR_UNBD')
    params.dualreductions = 0;
    result = gurobi(model, params);
    if strcmp(result.status, 'INFEASIBLE')
        exitflag = -2;
    elseif strcmp(result.status, 'UNBOUNDED')
        exitflag = -3;
    else
        exitflag = nan;
    end
elseif strcmp(result.status, 'INFEASIBLE')
    exitflag = -2;
elseif strcmp(result.status, 'UNBOUNDED')
    exitflag = -3;
else
    exitflag = nan;
end


if isfield(result, 'x')
    x = result.x;
else
    x = nan(n,1);
end

if isfield(result, 'objval')
    fval = result.objval;
else
    fval = nan;
end
