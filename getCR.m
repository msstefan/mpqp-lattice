function [ZA, zA, check] = getCR_v2(Q, Ht, c, A, b, F, active)

% The getCRandSol function uses the POP-QP representation to compute ...
%
% :params ...
% :returns: [ZA, zA]

Qi = inverse(Q);

if all(isnan(active))
%     SA = [];
%     sA = [];
%     LA = -0.5*Qi*Ht;
%     lA = -0.5*Qi*c;
    ZA = -0.5*A*Qi*Ht-F;
    zA = b+0.5*A*Qi*c;
    check = 1;
else
    inactive = 1:size(A,1);
    inactive(active(isnan(active) == 0)) = [];
    active(isnan(active))=[];

    KKTm = [Q A(active,:)'; A(active,:) zeros(numel(active))];

    if rank(KKTm) ~= size(KKTm,2)
        ZA = [];
        zA = [];
        check = 0;
        return
    end

    %
    Ai = inverse(A(active,:)*Qi * A(active,:).');
    %
    SA = -Ai * (2*F(active,:) + A(active,:)*Qi * Ht);
    sA = -Ai * (2*b(active,:) + A(active,:)*Qi * c);
    %
    LA = 0.5 * Qi*(A(active,:).' * Ai * (2*F(active,:) + A(active,:)*Qi * Ht) - Ht);
    lA = 0.5 *Qi*(A(active,:).' * Ai * (2*b(active) + A(active,:)*Qi * c) - c);
    %
    ZA = [A(inactive,:)*LA - F(inactive,:); -SA];
    zA = [b(inactive) - A(inactive,:)*lA; sA];
    %
    check = 1;
end