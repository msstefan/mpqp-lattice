function G_cal = minimal_sets(H,A,n)

index   = (1:n).';
V_H     = index;
V_H(H)  = [];
G_cal   = false(n);

V_H_label   = -ones(length(V_H),1); % -1: candidate, 0: min, 1: eliminate

G_column = 1;
for i = 1:length(V_H)
    v = V_H(i);
    h = closure_binary([H; v], A);

    if isempty(h)
        V_H_label(V_H == v) = 1;
        continue
    end

    h_v = setdiff(h, [H; v]);

    if ~isempty(h_v)
        m = length(h_v);

        flag = 0;

        for j = 1:m
            if V_H_label(V_H == h_v(j)) <= 0
                V_H_label(V_H == v) = 1;
                flag = 1;
                break
            end
        end

        if ~flag
            V_H_label(V_H == v) = 0;
            G_cal(h, G_column) = true;
            G_column = G_column + 1;
        end

    else
        V_H_label(V_H == v) = 0;
		
        G_cal(h, G_column) = true;
        G_column = G_column + 1;
    end
end

G_cal = unique(G_cal(:, any(G_cal,1)).', 'rows').';

end