% 相较 moveVertexOnMesh2 ，生成的网格是封闭网格了
% 和moveVertexOnMesh区别
% 1.在push across时，通过投影的方式，直接计算出了交点位置，以及跨越了那一条边
% 2.这样更加准确的实现算法

% 把每个点对应的面找到，最后统一做 delaunayTriangulation
function [vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh_new(vertices, ...
    vertices_cand, faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face, number_move)

% dam_ = 0.15 - 0.007*number_move; % k
nf = size(faces,1);
mesh.nf = nf;
nCand = length(nameF_cand);  
dam_ = 0.15;
% 计算新点的位置
vertices_cand_old = vertices_cand;
vertices_cand = vertices_cand + dam_ * forces; % P’ = P + kS

% 1.判断 push 不在原三角形内部的点
flag = isOnTriangle(nameF_cand, vertices_cand, vertices, faces, norm_face);

% 2.定义到网格三角形上
idx_out =  find(flag == 0);% 被push 到三角形外的点的顶点索引
num_out = length(idx_out);
if num_out  == 0
    return;
end
do_tri_idx = zeros(1,size(faces,1));
%% 1.寻找每个candidate 所在面和位置
for i = idx_out % i就是candidate 的索引
    push_point = vertices_cand_old(i,:);
    push_point_towards = vertices_cand(i,:);
    
    for k = 1:20
        v1 = vertices(faces(nameF_cand(i), 1),:);
        v2 = vertices(faces(nameF_cand(i), 2),:);
        v3 = vertices(faces(nameF_cand(i), 3),:);
        vs = [v1; v2; v3]; vs2 = [v2; v3; v1];
        [t1, t2] = solve_two_cross_lines(vs, vs2, ...
            repmat(push_point, 3, 1), repmat(push_point_towards, 3, 1));
        t1_flag = (t1>0 & t1<1); t2_flag = (t2>0 & t2<1);
        t_idx = find(t1_flag & t2_flag); % 找第一个与三角形相交点
        if length(t_idx) ~= 1 % 交点只能是0，所以停止 push across
            vertices_cand(i, :) = project_point_to_triangle(...
                push_point, vertices(faces(nameF_cand(i),:),:), norm_face(nameF_cand(i),:));
            break;
        end
        t2 = t2(t_idx);
        cross_point = push_point + t2*(push_point_towards - push_point);
        
        % 然后找旋转轴，更新 push_point
        e_num = t_idx;
        ek = [e_num, e_num+1];
        if e_num == 3
            ek = [3,1];
        end
        
        % 利用 ek 确定 push 到的下一个面
        idx_vk = faces(nameF_cand(i),[ek(2), ek(1)]);
        fnum_changed = hedge_face(idx_vk(1), idx_vk(2));
        
        % 绕 ek_inv 旋转 vout
        rotate_p = push_point_towards; % 旋转点
        rotate_line = vs(ek(2), :) - vs(ek(1), :);
        rotate_v0 = vs(ek(2), :);
        dot_ = dot(norm_face(nameF_cand(i),:), norm_face(fnum_changed,:));
        dot_(dot_ > 1) = 1; dot_(dot_ < -1) = -1;
        rotate_theta = acos(dot_);
        
        vp1 = ...
            point_rotate_line(rotate_p, rotate_line, rotate_v0, rotate_theta);
        vp2 = ...
            point_rotate_line(rotate_p, rotate_line, rotate_v0, - rotate_theta);
        vps = [vp1; vp2];
        
        dis_2 = dot(vps - rotate_v0, repmat(norm_face(fnum_changed,:),2,1), 2);
        dis_2 = abs(dis_2);
        idx_smaller = 1;
        if dis_2(1) > dis_2(2)
            idx_smaller = 2;
        end
        
        push_point_towards = vps(idx_smaller, :);
        push_point_old = push_point;
        push_point = cross_point + (1e-1)*(push_point_towards - cross_point);
        % push出去后，检查是否在 fnum_changed 所在的平面
        push_flag = isOnTriangle(fnum_changed, push_point, vertices, faces, norm_face);
        if ~push_flag
            vertices_cand(i, :) = project_point_to_triangle(...
                push_point_old, vertices(faces(nameF_cand(i),:),:), norm_face(nameF_cand(i),:));
            break;
        end
        vertices_cand(i, :) = push_point;
        do_tri_idx(nameF_cand(i)) = 1;
        nameF_cand(i) = fnum_changed;
        do_tri_idx(fnum_changed) = 1;
    end
% 所有candidate 所在面和位置得到更新
end

%% 2.根据面索引，将得到的网格重新参数化
% 找到没有添加和移除的三角形面上的连接关系
%每个面最多十个点，N = 10;
N = 10;
np = size(vertices,1); 
% 1.保持原来三角化，不重新做三角化的网格
faces_Mutual_old = zeros((2*N + 1)*size(faces,1), 3);
tm = 1;
for i = find(do_tri_idx == 0)
    idxs = find(nameF_cand == i) + np;
    num_idxs = length(idxs);
    % a.三角形内没有插入点
    if num_idxs == 0 
        faces_Mutual_old(tm,:) = faces(i,:);
        tm = tm + 1;
        continue;
    end
    % b.三角形内有一个插入点
    if num_idxs == 1
        [row, ~] = find(faces_Mutual == idxs);
        faces_Mutual_old(tm:(tm+2),:) = faces_Mutual(row,:);
        tm = tm+3;
        continue;
    end
    % c.三角形内有 2 ~ 10个插入点
    rows = zeros(1, (2*N + 1)*N); t = 1;
    for j = idxs
        [row, ~] = find(faces_Mutual == j);
        num_row = length(row);
        rows(t:(t+num_row-1)) = row;
        t = t + num_row;
    end
    rows = unique(rows(1:t-1));
    faces_rem = faces_Mutual(rows,:);
    num_rows = length(rows);
    faces_Mutual_old(tm:(tm+num_rows-1),:) = faces_rem;
    tm = tm + num_rows;
end
faces_Mutual_old = faces_Mutual_old(1:tm-1,:);

%2.改变了原来的三角化，重新做三角化的网格
faces_Mutual_new = zeros((2*N + 1)*size(faces,1), 3);
t = 1;
for i = find(do_tri_idx == 1)
    idx1 = faces(i,:);
    idx2 = find(nameF_cand == i);   
    vert1 = vertices(idx1,:); 
    vert2 = vertices_cand(idx2,:);
    idx = [idx1, idx2 + np];
    vert = [vert1; vert2];
    
    x_data = vert(:,1); y_data = vert(:,2);
    T = delaunayTriangulation(x_data, y_data);
    f_add = idx(T.ConnectivityList);
    % 保证生成的 mesh 是封闭定向一致的
    if norm_face(i,3) < 0
        f_add = f_add(:,[2,1,3]);
    end
    faces_Mutual_new(t:(t+size(f_add,1) - 1),:) = f_add;
    t = t + size(f_add,1);
end

faces_Mutual_new = faces_Mutual_new(1:t - 1,:);
faces_Mutual = [faces_Mutual_old; faces_Mutual_new];
end
