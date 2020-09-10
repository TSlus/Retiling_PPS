function [wights, tri] = uniform_point_in_triangle(N)

% clear; clc; N = 1;
points = zeros(N,2);
n = 0;
while (n < N)
    n = n + 1;
    x_i = 2 * rand - 1; % generate a number between -1 and 1
    y_i = 2 * rand; % generate a number between 0 and 2
    if (y_i >= (2 - 2 * abs(x_i)) - 0.1)  % if the points are outside the triangle
        n = n - 1; % decrease the counter to try to generate one more point
    elseif abs(y_i) <= 0.1
        n = n-1;
    else % if the point is inside the triangle
        points(n,:) = [x_i, y_i]; % add it to a list of points
    end
end
a = [0,0;1,0;0.5,1];
points(:,1) = points(:,1) + 1; points = points/2;

% figure(1)
% trimesh([1,2,3],a(:,1),a(:,2)); hold on;
% plot(points(:,1), points(:,2), 'r.');
% title ('6 points randomly distributed inside a triangle');

wights = zeros(N, 3);
v1 = repmat(a(3,:) - a(2,:), N, 1);
v2 = repmat(a(3,:) - a(1,:), N, 1);
v3 = repmat(a(2,:) - a(1,:), N, 1);

for1_cross = zeros(N, 3); for2_cross = zeros(N, 3); 
for1_cross(:,1:2) = v1; for2_cross(:,1:2) = points - a(2,:);  
cross1 = cross(for1_cross, for2_cross);
for1_cross(:,1:2) = points - a(1,:); for2_cross(:,1:2) = v2;  
cross2 = cross(for1_cross, for2_cross);
for1_cross(:,1:2) = v3; for2_cross(:,1:2) = points - a(1,:);  
cross3 = cross(for1_cross, for2_cross);

wights(:, 1) = sum(abs(cross1).^2, 2).^(1/2);
wights(:, 2) = sum(abs(cross2).^2, 2).^(1/2);
wights(:, 3) = sum(abs(cross3).^2, 2).^(1/2);

points = [a; points];
tri = delaunay(points(:,1),points(:,2));

% figure(2)
% trimesh(tri, points(:,1),points(:,2))

end