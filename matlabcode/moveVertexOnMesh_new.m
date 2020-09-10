% ��� moveVertexOnMesh2 �����ɵ������Ƿ��������
% ��moveVertexOnMesh����
% 1.��push acrossʱ��ͨ��ͶӰ�ķ�ʽ��ֱ�Ӽ�����˽���λ�ã��Լ���Խ����һ����
% 2.��������׼ȷ��ʵ���㷨

% ��ÿ�����Ӧ�����ҵ������ͳһ�� delaunayTriangulation
function [vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh_new(vertices, ...
    vertices_cand, faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face, number_move)

% dam_ = 0.15 - 0.007*number_move; % k
nf = size(faces,1);
mesh.nf = nf;
nCand = length(nameF_cand);  
dam_ = 0.15;
% �����µ��λ��
vertices_cand_old = vertices_cand;
vertices_cand = vertices_cand + dam_ * forces; % P�� = P + kS

% 1.�ж� push ����ԭ�������ڲ��ĵ�
flag = isOnTriangle(nameF_cand, vertices_cand, vertices, faces, norm_face);

% 2.���嵽������������
idx_out =  find(flag == 0);% ��push ����������ĵ�Ķ�������
num_out = length(idx_out);
if num_out  == 0
    return;
end
do_tri_idx = zeros(1,size(faces,1));
%% 1.Ѱ��ÿ��candidate �������λ��
for i = idx_out % i����candidate ������
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
        t_idx = find(t1_flag & t2_flag); % �ҵ�һ�����������ཻ��
        if length(t_idx) ~= 1 % ����ֻ����0������ֹͣ push across
            vertices_cand(i, :) = project_point_to_triangle(...
                push_point, vertices(faces(nameF_cand(i),:),:), norm_face(nameF_cand(i),:));
            break;
        end
        t2 = t2(t_idx);
        cross_point = push_point + t2*(push_point_towards - push_point);
        
        % Ȼ������ת�ᣬ���� push_point
        e_num = t_idx;
        ek = [e_num, e_num+1];
        if e_num == 3
            ek = [3,1];
        end
        
        % ���� ek ȷ�� push ������һ����
        idx_vk = faces(nameF_cand(i),[ek(2), ek(1)]);
        fnum_changed = hedge_face(idx_vk(1), idx_vk(2));
        
        % �� ek_inv ��ת vout
        rotate_p = push_point_towards; % ��ת��
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
        % push��ȥ�󣬼���Ƿ��� fnum_changed ���ڵ�ƽ��
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
% ����candidate �������λ�õõ�����
end

%% 2.���������������õ����������²�����
% �ҵ�û����Ӻ��Ƴ������������ϵ����ӹ�ϵ
%ÿ�������ʮ���㣬N = 10;
N = 10;
np = size(vertices,1); 
% 1.����ԭ�����ǻ��������������ǻ�������
faces_Mutual_old = zeros((2*N + 1)*size(faces,1), 3);
tm = 1;
for i = find(do_tri_idx == 0)
    idxs = find(nameF_cand == i) + np;
    num_idxs = length(idxs);
    % a.��������û�в����
    if num_idxs == 0 
        faces_Mutual_old(tm,:) = faces(i,:);
        tm = tm + 1;
        continue;
    end
    % b.����������һ�������
    if num_idxs == 1
        [row, ~] = find(faces_Mutual == idxs);
        faces_Mutual_old(tm:(tm+2),:) = faces_Mutual(row,:);
        tm = tm+3;
        continue;
    end
    % c.���������� 2 ~ 10�������
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

%2.�ı���ԭ�������ǻ������������ǻ�������
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
    % ��֤���ɵ� mesh �Ƿ�ն���һ�µ�
    if norm_face(i,3) < 0
        f_add = f_add(:,[2,1,3]);
    end
    faces_Mutual_new(t:(t+size(f_add,1) - 1),:) = f_add;
    t = t + size(f_add,1);
end

faces_Mutual_new = faces_Mutual_new(1:t - 1,:);
faces_Mutual = [faces_Mutual_old; faces_Mutual_new];
end
