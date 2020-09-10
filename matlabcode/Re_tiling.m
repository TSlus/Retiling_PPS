clear;clc;
load('model.mat');
% load('bronze.mat');

figure(1)
trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
title('原始网格');
%% 0.预计算
pre_compute; 

%% 1.加入Candidate vertices，得到Mutual tessellation
% nCand = nf; Candidate_vertices; 
nCand = nf; Candidate_vertices2; 

%% 
radius = 2 * sqrt(sa / nCand);
mesh.radius = radius;

figure(4) %到第37行，展示点的运动轨迹
plot3(vertices_cand(:,1),vertices_cand(:,2),vertices_cand(:,3),'b.');
vc_old = vertices_cand;
for number_move = 1:5
%% 2.计算每个Candidate vertices之间的距离，和排斥力，得到新的位置
forces = updata_candidate_vertices2(mesh, nameF_cand, vertices_cand);

nameF_cand2 = nameF_cand;
%% 3.通过每个点的受“力”，计算插入点新的位置
[vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh3(vertices, vertices_cand,...
    faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face);

end
force_test;

track_of_push;
% % 多次 push 后的网格 
% figure(5)
% trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));axis equal;
% title('push 后网格');

disp('=========================');
disp('push 结果测试：');
mesh_test;
disp('push 结果测试结束');
disp('=========================');

%% 4.Removing Old Vertices
% removing_old_vertices2_test; %相较于removing_old_vertices2除去了setdiff部分
removing_old_vertices3;
disp('=========================');
disp('Remove 结果测试：');
mesh_test;
disp('Remove 结果测试结束');
disp('=========================');
%% Remove 后图像
figure(6)
trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));axis equal;
title('Remove 后网格');
%% PPS流形--取最近点
faces_final = faces_Mutual;
vertices_final = vertices_Mutual;
load('surfaceP_PPS.mat');
figure(7)
plot3(surfaceP_PPS(:,1),surfaceP_PPS(:,2),surfaceP_PPS(:,3),'.');axis equal;
title('流形曲面');
for i = 1:nCand
    [~, idx] = min(sum(abs(vertices_final(i,:) - surfaceP_PPS).^2, 2));
    vertices_final(i,:) = surfaceP_PPS(idx, :);
end

%% 结合图像
figure(8)
trimesh(faces_final, vertices_final(:,1), vertices_final(:,2), vertices_final(:,3));axis equal;
title('Retiling & PPS');

%% retiling_bronze
% figure(9)
% load('retiling_bronze.mat');
% trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));axis equal;
% title('retiling bronze');