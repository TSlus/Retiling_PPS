function forces = updata_candidate_vertices2(mesh, nameF_cand, vertices_cand)
nCand = length(nameF_cand);    
forces = zeros(nCand, 3);

for i = 1:nCand
    % 判断的顶点是第 i 个面上顶点 i+ np
    % 1.直接计算同平面点的排斥力
    % 2.计算相邻面的排斥力，绕轴旋转
    % 3.较远的点，直接旋转到同平面
    
    nameF_cand_i = nameF_cand(i);
    % 执行1
    others1 = find(nameF_cand == nameF_cand_i);
    forces1 = vertices_cand(i, :) - vertices_cand(others1, :);
    
    % 执行2，与 i+np 有公共边的面
    vertices = mesh.vertices; faces = mesh.faces;
    hedge_face = mesh.hedge_face;
    norm_face = mesh.norm_face; radius = mesh.radius;
    
    a = faces(nameF_cand_i,1); b = faces(nameF_cand_i,2); c = faces(nameF_cand_i,3);
    f1 = hedge_face(b,a); f2 = hedge_face(c,b); f3 = hedge_face(a,c);
    % 待旋转的点
    others21 = find(nameF_cand == f1); n_1 = length(others21);
    others22 = find(nameF_cand == f2); n_2 = length(others22);
    others23 = find(nameF_cand == f3); n_3 = length(others23);
    others2 = [others21, others22, others23];
    
    p3 = vertices_cand(others2,:);
    % 旋转轴，直线上的点v0
    v0 = [vertices(a,:); vertices(b,:); vertices(c,:)];
    line = [vertices(b,:) - vertices(a,:); vertices(c,:) - vertices(b,:); vertices(a,:) - vertices(c,:)];
    % 旋转角度
    norm0 = repmat(norm_face(nameF_cand_i,:), 3, 1); % 第i个点所在面的法向量
    norms = [norm_face(f1,:); norm_face(f2,:); norm_face(f3,:)];
    dot_ = dot(norm0', norms');
    dot_(dot_ > 1) = 1; dot_(dot_ < -1) = -1;
    thetas = acos(dot_)';
    
    % 公式计算
    p3_rot_1 = point_rotate_line(p3, line, v0, - thetas, n_1, n_2, n_3); % p3旋转到 p 所在的平面
    p3_rot_2 = point_rotate_line(p3, line, v0,  thetas, n_1, n_2, n_3);
    
    norm_fi = norm_face(nameF_cand_i,:);
    p = vertices_cand(i,:);
    temp = [ones(1,n_1), 2*ones(1,n_2), 3*ones(1,n_3)];
    v0_detil = v0(temp, :);
    % 
    dist1 = dot(p3_rot_1 - v0_detil, repmat(norm_fi,(n_1 + n_2 + n_3),1), 2);
    dist2 = dot(p3_rot_2 - v0_detil, repmat(norm_fi,(n_1 + n_2 + n_3),1), 2);
    p3_cell{1} = p3_rot_1; p3_cell{2} = p3_rot_2;
    
    [~, min_idx] = min(abs([dist1, dist2]),[],2); %每一行的最小值和所在列
    p3_rot = p3_rot_1;
    for kk = 1:size(p3_rot,1)
        p3_rot(kk,:) = p3_cell{min_idx(kk)}(kk, :);
    end
    
    % 将p3_rot投影到三角形所在平面
    p3_rot = project_point_to_triangle(p3_rot, v0, norm_fi);
    % 计算到 p 的有向距离
    forces2 = p - p3_rot;
    
    % 执行3，先除去空间距离很远的点，剩下的点投影到平面上
    done_compute = zeros(1,nCand); done_compute([others1, others2]) = 1;
    others3 = find(done_compute == 0);
    dist_others3 = sum(abs(vertices_cand(others3,:) - p).^2, 2)';
    others3 = others3(dist_others3 < (1.5*radius)^2);
    
    forces3 = p - project_point_to_triangle(vertices_cand(others3,:), v0, norm_fi);
    
    force_p = [forces1; forces2; forces3];
    force_norm = sum(abs(force_p).^2, 2).^(1/2);
    force_idx = ((force_norm < radius) & (force_norm > 0));
    
    force_p = force_p(force_idx, :);
    force_norm = force_norm(force_idx, :);
    force_p = force_p./ force_norm ; % 将力单位化成--力方向
    force_p = sum((radius - force_norm) .* force_p, 1); %？？？ .*矩阵维度不相同
    forces(i,:) = force_p;
    
end

end
% 函数结束
