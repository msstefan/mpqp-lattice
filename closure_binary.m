function V = closure_binary(S, A)
%#codegen

V = [];
F = all(A(:,S), 2);

if isempty(F) || ~any(F == 1)
    return
end

V = all(A(F == 1,:), 1);
V = find(V == 1).';

end