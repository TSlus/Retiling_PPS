function [vertices_Mutual, faces_Mutual, nameF_cand, remain_p] = myRe_tiling(vertices, faces, n_k) %#ok<STOUT>
%% 0.Ԥ����
pre_compute_tiling;

%% 1.����Candidate vertices���õ�Mutual tessellation
% nCand = floor(nf/4); % �ײ��ӽ�����
nCand = floor(nf*n_k);
Candidate_vertices;

%%
radius = 2 * sqrt(sa / (nCand + size(vertices, 1)));
mesh.radius = radius;
for number_move = 1:5
    %% 2.����ÿ��Candidate vertices֮��ľ��룬���ų������õ��µ�λ��
%     forces = updata_candidate_vertices2(mesh, nameF_cand, vertices_cand);
    forces = updata_candidate_vertices_new(mesh, nameF_cand, vertices_cand);
    
    %% 3.ͨ��ÿ������ܡ����������������µ�λ��
%     [vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh3(vertices, vertices_cand,...
%         faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face, number_move);
      [vertices_cand, faces_Mutual, nameF_cand] = moveVertexOnMesh_new(vertices, vertices_cand,...
        faces, faces_Mutual, norm_face, forces, nameF_cand, hedge_face, number_move);

end

%% 4.Removing Old Vertices
removing_old_vertices2;
% removing_old_vertices_new;

end