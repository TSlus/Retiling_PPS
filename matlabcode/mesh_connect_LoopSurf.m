%% 原网格的LoopSurface（不同于细分两次后的）

% 1、先利用网格生成LoopSurface（原网格需要细分两次）
% 2、找到原网格三角形和细分两次后三角形之间的关系（hedge_face）
% 3、得到原网格三角形和LoopSurface点关系（原网格的LoopSurface）
function loop_point = mesh_connect_LoopSurf(vertices, faces)
%% 1.LoopSurfaces
% load('model.mat');
faces_old = faces;
% 先做两次细分
loopTimes = 2;
[vertices, faces, hedge_face] = myLoop(vertices, faces, loopTimes);
% trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)); axis equal

% 生成loopSurface
% 注：这里的loop_point与细分后三角形对应
loop_point = LoopSurfaceSample(vertices, faces); 

%% 2.通过原网格三角形和细分两次后三角形之间的关系
% 得到原网格的LoopSurface

nf = size(faces_old, 1);

x1 = faces(:,1); x2 = faces(:,2); x3 = faces(:,3);
X = [x1; x2; x3]; Y = [x2; x3; x1];
R = [1:size(faces, 1), 1:size(faces, 1), 1:size(faces, 1)]';
he_f = sparse(X, Y, R);
connect{nf} = []; % nf 是细分前faces个数
nEvery = 100; % loop细分曲面loop_point每个三角形采点最多100个
for i=1:nf
    [row, col] = find(hedge_face == i);
    fi = he_f(sub2ind(size(he_f), row, col));
    fi = unique(fi); n_fi = length(fi);
    cell_p = loop_point(fi);
    % 将这些元胞中的点取出来
    v_fi = zeros(n_fi*nEvery, 3);
    t = 1;
    for j = 1:n_fi
        k = size(cell_p{j}, 1);
        v_fi(t:t+k-1, :) = cell_p{j};
        t = t+k;
    end
     v_fi =  v_fi(1:t-1, :);
     connect{i} = v_fi;
end
% 注：这里的loop_point与细分前三角形对应
loop_point = connect; 

end
