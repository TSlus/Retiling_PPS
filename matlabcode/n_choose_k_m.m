% n_choose_k
% n最大值为(mu + 1)

clear;clc;
% load('sbunnySimp.mat');
load('bronze.mat');

numP = size(vertices,1);numF = size(faces,1);
% 计算顶点的邻面adjFaces索引
v_adjFaces{numP} = [];
for i = 1 : numP
    [v_adjFaces{i},~] = find(faces == i) ;%第几个面与定点i相关
end

% 每个顶点的度v_valence
v_valence(numP)=0;
for i=1:numP
    [v_valence(i),~] = size(v_adjFaces{i});
end
max_val = max(v_valence);
n_choose_k = sparse(max_val + 2, max_val + 2);

for i = 0:max_val + 1
    for j = 0:max_val + 1
    if i >= j
        n_choose_k(i+1, j+1) = nchoosek(i,j);
    end
    end
end