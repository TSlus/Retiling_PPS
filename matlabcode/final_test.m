clear;clc
load('model.mat');
%% 1.���� loopSurface 
loop_point = mesh_connect_LoopSurf(vertices, faces);

%% 2.PPS
surfaceP = surfaceConstruction(vertices, faces, loop_point);
