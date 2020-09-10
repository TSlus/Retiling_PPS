% 将移动 candidate point 后得到的 faces_Mutual，去除 old_point
vertices_Mutual = [vertices; vertices_cand];
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
        
    hedge(i, :) = 0;
    hedge(:, i) = 0;
    
    x1 = faces_add2(:,1); x2 = faces_add2(:,2); x3 = faces_add2(:,3);
    X = [x1; x2; x3]; Y = [x2; x3; x1];
    hedge(sub2ind(size(hedge), X, Y)) = 1;
    nhe = size(X, 1);
    if (mod(nhe_old - 3*num_nearp + nhe, 2) ~= 0)
        warning('(mod(nhe, 2) ~= 0)');
    end
    
    if (full(sum(sum(hedge))) ~= (nhe_old - 3*num_nearp + nhe))
        faces_Mutual = faces_Mutual_old;
        hedge = hedge_old;
        continue;
    end
    nhe_old = nhe_old - 3*num_nearp + nhe;
    
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

remove_mat = sparse(sparse_idx(:,1),sparse_idx(:,2),-1);
idx_adjust = full(sum(remove_mat, 1));
idx_adjust = (1:np) + idx_adjust;
idx_adjust = [idx_adjust, ((np+1):(np+nCand)) - length(remove_p)];
faces_Mutual = idx_adjust(faces_Mutual);
