% �ҵ� old_point ��Ӧ�������������


% ����ÿ�� candidate_point ��Ȩֵ
v1 = vertices(faces(nameF_cand, 1),:);
v2 = vertices(faces(nameF_cand, 2),:);
v3 = vertices(faces(nameF_cand, 3),:);

norm_fi = norm_face(nameF_cand, :);
a1 = cross(v3 - v2, vertices_cand - v2, 2);
a2 = cross(vertices_cand - v1, v3 - v1, 2);
a3 = cross(v2 - v1, vertices_cand - v1, 2);

dot1 = dot(a1', norm_fi',1); dot2 = dot(a2', norm_fi', 1); dot3 = dot(a3', norm_fi',1);
cand_flag = (dot1 > 0) & (dot2 > 0) & (dot3 > 0);

if isempty(find(cand_flag == 0,1))
    disp('ÿ��candidate point ��nameF_cand �ڲ�!')
else
    disp('����candidate point ����nameF_cand �ڲ�.');
    return;
end

% alpha = sum(abs(a1).^2, 2).^(1/2) ./ (2*si); 
% beta = sum(abs(a2).^2, 2).^(1/2) ./ (2*si); 
% gamma = sum(abs(a3).^2, 2).^(1/2) ./ (2*si);

% ����Ȩ��ֱ���������������ȥ��

