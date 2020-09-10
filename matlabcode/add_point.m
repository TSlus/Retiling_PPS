% 第一条半边情况
if flag_he(i) % 如果第一条半边上已经插入了点
    k1 = flag_he(i);
else          % 如果第一条半边上还没插入新点
    k1 = idxp;
    flag_he(i) = idxp;
    % 对面半边也要插入相应的点
    name = hedge_oppo_name(a, b);
    flag_he(name) = idxp;
    idxp = idxp + 1;
end

% 第二条半边情况
if flag_he(i + nf) % 如果第二条半边上已经插入了点
    k2 = flag_he(i + nf);
else          % 如果第二条半边上还没插入新点
    k2 = idxp;
    flag_he(i + nf) = idxp;
    % 对面半边也要插入相应的点
    name = hedge_oppo_name(b, c);
    flag_he(name) = idxp;
    idxp = idxp + 1;
end

% 第三条半边情况
if flag_he(i + 2*nf) % 如果第一条半边上已经插入了点
    k3 = flag_he(i + 2*nf);
else          % 如果第一条半边上还没插入新点
    k3 = idxp;
    flag_he(i + 2*nf) = idxp;
    % 对面半边也要插入相应的点
    name = hedge_oppo_name(c, a);
    flag_he(name) = idxp;
    idxp = idxp + 1;
end