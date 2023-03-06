clc, clear, close all

% data = load("mpQP_test_2.mat");
data = load("example_data.mat");

% Create the lifted domain
ldom = Polyhedron([data.problem.A -data.problem.F], data.problem.b);

tic
V_polar = (ldom.A ./ ldom.b).';
P_polar = Polyhedron(V_polar.');

% Compute the mpQP solution using the modified face lattice algorithm
tol = 1e-6;

% vertex-facet incidence matrix
A   = (abs(P_polar.A * V_polar - P_polar.b) <= tol);

% initialize the 'queue' and the solution stack L
Q{1} = [];
L = digraph;
L = addnode(L, cellstr(strjoin(string([]))));
n = size(V_polar,2);
solutions = 1;

while ~isempty(Q)
    % Extract the 1st element from Q
    H = Q{1};
    % Remove the element from Q
    Q(1) = [];
    % Compute minimal sets
    G = minimal_sets(H,A,n);
    % Compute size
    G_len = size(G,2);
    % Find a solution
    for i = 1:G_len
        S = find(G(:,i) == 1);
        closure = closure_binary(S,A);
        element = cellstr(strjoin(string(closure)));
        % Search the digraph
        found = findnode(L, element);

        if ~found
            [ZA, zA, check] = getCR(data.problem.Q, data.problem.Ht, data.problem.c, data.problem.A, data.problem.b, data.problem.F, closure);
            % Test the rank condition
            if ~check
                continue
            end
            % Check if empty
            if Polyhedron(ZA,zA).isEmptySet
                continue
            end
            % If the rank condition passes and the CR is not empty then
            % add to queue
            Q{end + 1} = closure;
            solutions = solutions + 1;
        end
        % add to solution graph
        L = addedge(L, cellstr(strjoin(string(H))), element);
    end
end
time = toc;

fprintf("Found %d solutions in %0.6f seconds!\n", solutions, time)

plot(L)

% ToDo: Solve with yalmip

