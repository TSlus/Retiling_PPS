function [vertices, faces, hedge_face] = LoopSubdivide(vertices, faces, hedge_face_old)
nf = size(faces, 1); np = size(vertices, 1);
nhe = 3 * nf; % �����
ne = nhe / 2;
if mod(nhe, 2) ~= 0
    warning('����������Ŀ����ż���������Ƿ������')
    return;
end

% ��߽ṹ
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; 
% ���
he = zeros(nhe, 2); 
he(:,1) = X; he(:,2) = Y;
% ���--��
F = [1:nf, 1:nf, 1:nf]';
hedge_face = sparse(X, Y, F);
% ���--�Ե�
Z = [x3; x1; x2];
hedge_third_p = sparse(X, Y, Z);
% ���--�����ߵ�����
hedge_oppo_name = sparse(Y, X, (1:3*nf)');

% ÿ��ԭʼ����Ķ�
v_valence(np) = 0;
nearPs{np} = [];
for P = 1:np
    neighbor_P = find(hedge_face(:,P));
    nearPs{P} = neighbor_P;           % P�������
    v_valence(P) = length(neighbor_P);  % P�Ķ�
end

% ˼·�����ѭ������ÿ������ѳ��ĸ�
flag_he = zeros(nhe, 1); % ��¼������Ƿ�������µ㣬ͬʱ������¼����ĵ�
% hedge_midpoint = hedge_face; % ��ÿ������ϲ���һ����

% ��ѭ��
idxp = np + 1; % ������µ�����
faces_new = zeros(4*nf, 3);
vertices_new = zeros(np + ne, 3);
flag_v = vertices_new(:, 1); % ��¼�����Ƿ��Ѿ�����
face_name = zeros(12*nf, 3); % ��¼�²�����ߣ�����ԭʼ��
for i = 1:nf
    a = faces(i, 1); b = faces(i, 2); c = faces(i, 3);
    % �����µ�k1, k2, k3
    add_point; 
    % 1.������
    face_add  = [a, k1, k3; 
                 k1, b, k2;
                  k3, k2, c;
                  k1, k2, k3];
    faces_new(4*i - 3:4*i,:) = face_add;
    
    y1 = [face_add(:,1),face_add(:,2)]; y2 = [face_add(:,2),face_add(:,3)]; 
    y3 = [face_add(:,3),face_add(:,1)]; 
    face_name((12*i-11):12*i, [1,2]) = [y1; y2; y3];
    face_name((12*i-11):12*i, 3) = hedge_face_old(a, b);
    
    % 2.���µ�
    update_vertex;
end

vertices = vertices_new;
faces = faces_new;
hedge_face = sparse(face_name(:, 1), face_name(:, 2), face_name(:, 3));
end


