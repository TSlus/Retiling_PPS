% do_triangulation;
% triangulation½øÐÐ
vertex_tri1 = vertices(vertex_tri_idx1, :);
vertex_tri2 = vertices_cand(vertex_tri_in, :);

vertex_tri_idx = [vertex_tri_idx1, vertex_tri_in + size(vertices,1)];
vertex_tri = [vertex_tri1; vertex_tri2];
x_data = vertex_tri(:,1); y_data = vertex_tri(:,2);
T = delaunayTriangulation(x_data, y_data); 
