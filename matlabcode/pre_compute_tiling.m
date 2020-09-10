nf = size(faces,1); np = size(vertices, 1);

%% 1.判断网格封闭性
x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1]; % R = [1:size(faces,1), 1:size(faces,1), 1:size(faces,1)]';
halfEdge = sparse(X, Y, 1);
if ~isempty(find(halfEdge - halfEdge', 1))
    warning('网格不是封闭网格。');
    return;
end

%% 2.利用半边来索引面
hedge_face = sparse(X, Y, [1:nf, 1:nf, 1:nf]');

%% 3.计算每个面的法向量
v1 = vertices(faces(:,2),:) - vertices(faces(:,1),:);
v2 = vertices(faces(:,3),:) - vertices(faces(:,1),:);
norm_face = cross(v1, v2, 2);
norm_face_normal = sum(abs(norm_face).^2, 2).^(1/2);
norm_face = norm_face ./ norm_face_normal;
% onee = sum(norm_face.^2,2) - 1;max(abs(onee))

%% 4.area of suface
si = norm_face_normal / 2; %每个小三角形面积
sa = sum(si);              %网格表面积
% 记录si过小的面
% 百分之一的“小面积”三角形上不选 candidate_points
k_si = ceil(nf * 0.05);
[~, si_small] = mink(si, k_si);

%% 将网格信息整体化
mesh.vertices = vertices; mesh.faces = faces;
mesh.np = np; mesh.nf = nf;
mesh.hedge_face = hedge_face;
mesh.norm_face = norm_face; 
mesh.si = si; mesh.sa = sa;

