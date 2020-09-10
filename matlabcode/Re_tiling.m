clear;clc;
load('model.mat');
% load('bronze.mat');

figure(1)
trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
title('ԭʼ����');
%% 0.Ԥ����
pre_compute; 

%% 1.����Candidate vertices���õ�Mutual tessellation
% nCand = nf; Candidate_vertices; 
nCand = nf; Candidate_vertices2; 

%% 
radius = 2 * sqrt(sa / nCand);
mesh.radius = radius;

figure(4) %����37�У�չʾ����˶��켣
plot3(vertices_cand(:,1),vertices_cand(:,2),vertices_cand(:,3),'b.');
vc_old = vertices_cand;
for number_move = 1:5
%% 2.����ÿ��Candidate vertices֮��ľ��룬���ų������õ��µ�λ��
forces = updata_candidate_vertices2(mesh, nameF_cand, vertices_cand);

nameF_cand2 = nameF_cand;
%% 3.ͨ��ÿ������ܡ����������������µ�λ��
[vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh3(vertices, vertices_cand,...
    faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face);

end
force_test;

track_of_push;
% % ��� push ������� 
% figure(5)
% trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));axis equal;
% title('push ������');

disp('=========================');
disp('push ������ԣ�');
mesh_test;
disp('push ������Խ���');
disp('=========================');

%% 4.Removing Old Vertices
% removing_old_vertices2_test; %�����removing_old_vertices2��ȥ��setdiff����
removing_old_vertices3;
disp('=========================');
disp('Remove ������ԣ�');
mesh_test;
disp('Remove ������Խ���');
disp('=========================');
%% Remove ��ͼ��
figure(6)
trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));axis equal;
title('Remove ������');
%% PPS����--ȡ�����
faces_final = faces_Mutual;
vertices_final = vertices_Mutual;
load('surfaceP_PPS.mat');
figure(7)
plot3(surfaceP_PPS(:,1),surfaceP_PPS(:,2),surfaceP_PPS(:,3),'.');axis equal;
title('��������');
for i = 1:nCand
    [~, idx] = min(sum(abs(vertices_final(i,:) - surfaceP_PPS).^2, 2));
    vertices_final(i,:) = surfaceP_PPS(idx, :);
end

%% ���ͼ��
figure(8)
trimesh(faces_final, vertices_final(:,1), vertices_final(:,2), vertices_final(:,3));axis equal;
title('Retiling & PPS');

%% retiling_bronze
% figure(9)
% load('retiling_bronze.mat');
% trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));axis equal;
% title('retiling bronze');