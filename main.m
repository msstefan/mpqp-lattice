clc, clear, close all

import java.util.LinkedList

% data = load("mpQP_test_2.mat");
data = load("example_data.mat");

% Create the queue
Q = LinkedList();

tic
% Create the lifted domain
V_polar = ([data.problem.A -data.problem.F] ./ data.problem.b).';

% Compute the mpQP solution using the modified face lattice algorithm
tol = 1e-6;

% vertex-facet incidence matrix
hull = convhulln(V_polar', {'Qt','Qx'});
nineq = size(hull,1);
A = false(nineq, size(V_polar,2));
% create a matrix of repeated column indices
colIndexMatrix = repmat((1:nineq)', 1, size(hull, 2));
% convert subscripts to linear indices and set the corresponding elements to true
A(sub2ind(size(A), colIndexMatrix, hull)) = true;

% initialize the 'queue' and the solution stack L
% Q{1} = [];
Q.add([]);
L = digraph;
L = addnode(L, cellstr(strjoin(string([]))));
n = size(V_polar,2);
solutions = 1;

while Q.size() > 0
    % Extract the 1st element from Q
    H = Q.remove();
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
            Q.add(closure);
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

