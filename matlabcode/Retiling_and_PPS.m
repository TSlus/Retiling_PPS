clear;clc;
%% 输入
% 1.网格
% load('model.mat');
load('bunny.mat');
% load('bronze.mat'); 

% 2.顶点数
n_k = 1/4; % 报告中 n_k 取1/4、1/2、1 * numF，作为选取的点数目

%%
figure(1)
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3));axis equal
title('初始网格')
%% do Re_Tiling
% retiling 最后要用new中加速

[vertices_Mutual, faces_Mutual, nameF_cand, remain_p] = myRe_tiling(vertices, faces, n_k);

figure(2)
trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));axis equal;
title('mesh after removing old points');

%%
% return; 
% 计算出每个点的重心坐标，再做PPS
% 减少流形计算部分的时间

%% do PPS
% 1.loopSurface
loop_point = mesh_connect_LoopSurf(vertices, faces); % #ok<UNRCH>
% PPS
surfaceP = surfaceConstruction(vertices, faces, loop_point);
figure(3)
plot3(surfaceP(:,1),surfaceP(:,2),surfaceP(:,3),'b.');axis equal
title('PPS')
%% Retiling and PPS
vertices_final = vertices_Mutual;faces_final = faces_Mutual;
for i = 1:size(vertices_Mutual)
    [~, idx] = min(sum(abs(vertices_final(i,:) - surfaceP).^2, 2));
    vertices_final(i,:) = surfaceP(idx, :);
end

%图像
figure(4)
trimesh(faces_final, vertices_final(:,1), vertices_final(:,2), vertices_final(:,3));axis equal;
title('Retiling & PPS');