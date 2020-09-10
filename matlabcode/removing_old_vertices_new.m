% 将移动 candidate point 后得到的 faces_Mutual，去除 old_point
remove_idx = zeros(1, np); % 被移除就赋值为1
nCand = size(vertices_cand, 1);
sparse_idx = zeros(ceil(np*np), 2); ts =1;

% 循环中修改半边结构，速度或许更快
edges =[faces_Mutual(:,[1,2]); faces_Mutual(:,[2,3]); faces_Mutual(:,[3,1])];
% 半边数
nhe_old = size(edges, 1);
x = edges(:,1); y = edges(:,2); z = ones(nhe_old, 1);
hedge = sparse(x, y, z);

for i = 1:np
    faces_Mutual_old = faces_Mutual;
    hedge_old = hedge;
    
    [r_i, ~] = find(faces_Mutual == i);
    delete_fM = zeros(1, size(faces_Mutual, 1));
    delete_fM(r_i) = 1;
    faces_i = faces_Mutual(r_i,:);
    near_p = findNearP(faces_i, i);
    num_nearp = length(near_p);
    
    t = (1:num_nearp) * 2 * pi / num_nearp; t = fliplr(t);
    x = cos(t); y = sin(t);
    xy = [x; y]; xy = xy';
    
    flip_near = fliplr(near_p); % 计算OneRing，顺时针给出
    di_vec = vertices_Mutual(flip_near,:) - vertices_Mutual(i,:);
    di = sum(abs(di_vec).^2, 2).^(1/2);
    
    xy = xy .* di;
    polyin = polyshape(xy(:,1), xy(:,2));
    T = triangulation(polyin);
    
    % 加入的面
    faces_add2 = flip_near(T.ConnectivityList);
    % 新增的面
    faces_add1 = faces_Mutual(delete_fM == 0, :);
    faces_Mutual = [faces_add1; faces_add2];
    
    % 把one-ring的半边去掉
    hE = zeros(num_nearp, 2);
    for j = 1:num_nearp - 1
        hE(j, :) = near_p(j:j+1); 
    end
    hE(end, :) = [near_p(end), near_p(1)];
    hedge(sub2ind(size(hedge), hE(:,1), hE(:,2))) = 0;
    
    x1 = faces_add2(:,1); x2 = faces_add2(:,2); x3 = faces_add2(:,3);
    X = [x1; x2; x3]; Y = [x2; x3; x1];
    hedge_Z = hedge(sub2ind(size(hedge), X, Y));
    hedge(sub2ind(size(hedge), X, Y)) = hedge_Z + 1;

    if max(max(hedge)) > 1
        faces_Mutual = faces_Mutual_old;
        hedge = hedge_old;
        continue;
    end
    
    remove_idx(i) = 1;
    sparse_idx(ts:(ts+np-i-1),1) = i;
    sparse_idx(ts:(ts+np-i-1),2) = ((i+1):np)';
    ts = ts + np-i;
    %     remove_mat(i, i+1:end) = -1; % 利用稀疏矩阵加速
end
sparse_idx = sparse_idx(1:ts-1,:);
remain_p = find(remove_idx == 0);
remove_p = find(remove_idx == 1);

vertices_Mutual = [vertices_Mutual(remain_p, :); vertices_cand];

% 为什么稀疏矩阵比是矩阵赋值更慢？？？？？？
% remove_mat = sparse(sparse_idx(:,1),sparse_idx(:,2),-1);
ia = max(sparse_idx(:,1)); ib = max(sparse_idx(:,2));
remove_mat = zeros(ia, ib);
remove_mat(sub2ind(size(remove_mat), sparse_idx(:,1),sparse_idx(:,2))) = -1;

idx_adjust = full(sum(remove_mat, 1));
idx_adjust = (1:np) + idx_adjust;
idx_adjust = [idx_adjust, ((np+1):(np+nCand)) - length(remove_p)];
faces_Mutual = idx_adjust(faces_Mutual);

faces_final = faces_Mutual;
vertices_final = vertices_Mutual;
