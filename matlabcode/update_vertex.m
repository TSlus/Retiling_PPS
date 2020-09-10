% aµã
if flag_v(a) == 0
    flag_v(a) = 1;
    va = vertices(a,:);
    u = v_valence(a);
    beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
    vertices_new(a,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{a},:);
end

% bµã
if flag_v(b) == 0
    flag_v(b) = 1;
    va = vertices(b,:);
    u = v_valence(b);
    beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
    vertices_new(b,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{b},:);
end

% cµã
if flag_v(c) == 0
    flag_v(c) = 1;
    va = vertices(c,:);
    u = v_valence(c);
    beta = 1/u*(5/8-(3/8+1/4*cos(2*pi/u))^2);
    vertices_new(c,:) = (1-u*beta)*va + ones(1,u) * beta * vertices(nearPs{c},:);
end

%%
% k1
if flag_v(k1) == 0
    flag_v(k1) = 1;
    v1 = vertices(a,:); v2 = vertices(b, :); v3 = vertices(c, :);
    d = hedge_third_p(b, a); v4 = vertices(d,:);
    vertices_new(k1, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
end

% k2
if flag_v(k2) == 0
    flag_v(k2) = 1;
    v1 = vertices(b,:); v2 = vertices(c, :); v3 = vertices(a, :);
    d = hedge_third_p(c, b); v4 = vertices(d,:);
    vertices_new(k2, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
end

% k3
if flag_v(k3) == 0
    flag_v(k3) = 1;
    v1 = vertices(c,:); v2 = vertices(a, :); v3 = vertices(b, :);
    d = hedge_third_p(a, c); v4 = vertices(d,:);
    vertices_new(k3, :) = (v1 + v2)*3/8 + (v3 + v4)/8;
end

