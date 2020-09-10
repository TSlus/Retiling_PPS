function forces = updata_candidate_vertices_new(mesh, nameF_cand, vertices_cand)
nCand = length(nameF_cand);   
forces = zeros(nCand, 3);
vertices = mesh.vertices; faces = mesh.faces;
nf = mesh.nf;
hedge_face = mesh.hedge_face;
norm_face = mesh.norm_face; radius = mesh.radius;

% ������Ҫ��ѭ��д��һ����

% ����һ��Ԫ�������ڸ�������Щ�㣬����find����ʹ��
% t: ÿһ����ָʾ������������������ڴ�
t = ones(1, mesh.nf);
tem = zeros(1, 10); % ÿ���棬������10��candidate��
cell_p(1:mesh.nf) = {tem}; 
for i = 1:nCand
    fi = nameF_cand(i);
    ti = t(fi); % ��fi�����ָ��
    cell_p{fi}(ti) = i;
    t(fi) = ti + 1;
end
% ȡ�����������
for i = 1:nf
    ti = t(i);
    cell_p{i} = cell_p{i}(1:(ti-1));
end

% 1.����ÿ�� candidate ��� forces1

for i = 1:nCand
    nameF_cand_i = nameF_cand(i);
    % ִ��1
    others1 = cell_p{nameF_cand_i};
    forces1 = vertices_cand(i, :) - vertices_cand(others1, :);
    
    % ִ��2���� i+np �й����ߵ���
    a = faces(nameF_cand_i,1); b = faces(nameF_cand_i,2); c = faces(nameF_cand_i,3);
    f1 = hedge_face(b,a); f2 = hedge_face(c,b); f3 = hedge_face(a,c);
    % ����ת�ĵ�
    others21 = cell_p{f1}; n_1 = length(others21);
    others22 = cell_p{f2}; n_2 = length(others22);
    others23 = cell_p{f3}; n_3 = length(others23);
    others2 = [others21, others22, others23];
    
    p3 = vertices_cand(others2,:);
    % ��ת�ᣬֱ���ϵĵ�v0
    v0 = [vertices(a,:); vertices(b,:); vertices(c,:)];
    line = [vertices(b,:) - vertices(a,:); vertices(c,:) - vertices(b,:); vertices(a,:) - vertices(c,:)];
    % ��ת�Ƕ�
    norm0 = repmat(norm_face(nameF_cand_i,:), 3, 1); % ��i����������ķ�����
    norms = [norm_face(f1,:); norm_face(f2,:); norm_face(f3,:)];
    dot_ = dot(norm0', norms');
    dot_(dot_ > 1) = 1; dot_(dot_ < -1) = -1;
    thetas = acos(dot_)';
    
    % ��ʽ����
    p3_rot_1 = point_rotate_line(p3, line, v0, - thetas, n_1, n_2, n_3); % p3��ת�� p ���ڵ�ƽ��
    p3_rot_2 = point_rotate_line(p3, line, v0,  thetas, n_1, n_2, n_3);
    
    norm_fi = norm_face(nameF_cand_i,:);
    p = vertices_cand(i,:);
    temp = [ones(1,n_1), 2*ones(1,n_2), 3*ones(1,n_3)];
    v0_detil = v0(temp, :);
    %
    dist1 = dot(p3_rot_1 - v0_detil, repmat(norm_fi,(n_1 + n_2 + n_3),1), 2);
    dist2 = dot(p3_rot_2 - v0_detil, repmat(norm_fi,(n_1 + n_2 + n_3),1), 2);
    p3_cell{1} = p3_rot_1; p3_cell{2} = p3_rot_2;
    
    [~, min_idx] = min(abs([dist1, dist2]),[],2); %ÿһ�е���Сֵ��������
    p3_rot = p3_rot_1;
    for kk = 1:size(p3_rot,1)
        p3_rot(kk,:) = p3_cell{min_idx(kk)}(kk, :);
    end
    
    % ��p3_rotͶӰ������������ƽ��
    p3_rot = project_point_to_triangle(p3_rot, v0, norm_fi);
    % ���㵽 p ���������
    forces2 = p - p3_rot;
    
    % ִ��3���ȳ�ȥ�ռ�����Զ�ĵ㣬ʣ�µĵ�ͶӰ��ƽ����
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
    force_p = force_p./ force_norm ; % ������λ����--������
    force_p = sum((radius - force_norm) .* force_p, 1); %������ .*����ά�Ȳ���ͬ
    forces(i,:) = force_p; 
end

end