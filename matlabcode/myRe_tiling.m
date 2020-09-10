function [vertices_Mutual, faces_Mutual, nameF_cand, remain_p] = myRe_tiling(vertices, faces, n_k) %#ok<STOUT>
%% 0.预计算
pre_compute_tiling;

%% 1.加入Candidate vertices，得到Mutual tessellation
% nCand = floor(nf/4); % 底部接近均匀
nCand = floor(nf*n_k);
Candidate_vertices;

%%
radius = 2 * sqrt(sa / (nCand + size(vertices, 1)));
mesh.radius = radius;
for number_move = 1:5
    %% 2.计算每个Candidate vertices之间的距离，和排斥力，得到新的位置
%     forces = updata_candidate_vertices2(mesh, nameF_cand, vertices_cand);
    forces = updata_candidate_vertices_new(mesh, nameF_cand, vertices_cand);
    
    %% 3.通过每个点的受“力”，计算插入点新的位置
%     [vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh3(vertices, vertices_cand,...
%         faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face, number_move);
      [vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh_new(vertices, vertices_cand,...
        faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face, number_move);

end

%% 4.Removing Old Vertices
removing_old_vertices2;
% removing_old_vertices_new;

end