% 计算每个顶点的绝对曲率
FV.vertices = vertices;
FV.faces = faces;
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2] = GetCurvature(FV, true);
abs_Cmean = abs(Cmean);
abs_Cmean_faces = (abs_Cmean(faces(:,1)) + abs_Cmean(faces(:,2)) + abs_Cmean(faces(:,3)));

k1 = ceil(nf * 0.2);
k3 = ceil(nf * 0.2);
k2 = nf - k1 - k3;
[~, old_face_name1] = mink(abs_Cmean_faces,k1);
[~, old_face_name3] = maxk(abs_Cmean_faces,k3);
old_face_name2 = setdiff(1:nf, [old_face_name1; old_face_name3])';

nEveryFace = (nCand - k2) / k3;
nEveryFace = ceil(nEveryFace); % 这样生成出来的candidate point 会比输入的多
if nCand < k3
    disp('the quantity of candidate points is too little.')
    return;
end
if nCand < k2 + k3
    old_face_name1 = [old_face_name1; old_face_name2(nCand - k3 + 1:end)];
    k1 = length(old_face_name1);
    old_face_name2 = old_face_name2(1:nCand - k3);  
    k2 = length(old_face_name2);
    nEveryFace = 1;
end

% 调整 old_face_name 中元素
adjust_cand_vertices;

tem = repmat(1:k3, nEveryFace, 1);
nameF_cand = [old_face_name2; old_face_name3(tem(:))]'; % 精彩！！！
nCand = length(nameF_cand);

% 执行步骤1
faces_Mutual_k1 = faces(old_face_name1, :);

% 执行步骤2，产生一个点
v_1 = vertices(faces(old_face_name2,1),:);
v_2 = vertices(faces(old_face_name2,2),:);
v_3 = vertices(faces(old_face_name2,3),:);
vertices_cand_k2 = 1/3 * (v_1 + v_2 + v_3);
add_vertex_name = np + k2; %加入点的名字
% 加入面
faces_Mutual_k2 = zeros(3*nf, 3);
ti = 1; 
for i = old_face_name2'
    faces_Mutual_k2(3*ti-2:(3*ti),:) = [faces(i,1), faces(i,2), ti+np;
                                 faces(i,2), faces(i,3), ti+np;
                                 faces(i,3), faces(i,1), ti+np  ];
    ti = ti + 1;
end
faces_Mutual_k2 = faces_Mutual_k2(1:3*(ti - 1),:);

% 执行步骤3
[wights, tri] = uniform_point_in_triangle(nEveryFace);
num_add_f = size(tri, 1);

vertices_cand_k3 = zeros(nf * nEveryFace, 3);
v_ti = 1;
faces_Mutual_k3 = zeros(nf * num_add_f, 3); % 每个三角形面剖分 num_add_f 个三角形
ti = 1;
for i = old_face_name3'
    v3_f = vertices(faces(i,:),:);
    vcand_new = wights * v3_f;
    vertices_cand_k3(v_ti:(v_ti + nEveryFace - 1), :) = vcand_new;
    v_ti = v_ti + nEveryFace;
    
    face_order = [faces(i,:), (add_vertex_name + 1):(add_vertex_name + nEveryFace)];
    add_vertex_name = add_vertex_name + nEveryFace;
    
    face_add = face_order(tri);
    faces_Mutual_k3(ti:(ti + num_add_f - 1),:) = face_add;
    ti = ti + num_add_f;
end
vertices_cand_k3 = vertices_cand_k3(1:v_ti - 1, :);
faces_Mutual_k3 = faces_Mutual_k3(1:ti - 1,:);

vertices_cand = [vertices_cand_k2; vertices_cand_k3];
faces_Mutual = [faces_Mutual_k1; faces_Mutual_k2; faces_Mutual_k3];

vertices_Mutual = [vertices; vertices_cand];




