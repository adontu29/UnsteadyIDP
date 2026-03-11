clear;close all; clc;
S    = load('U15A10-10K100P0000.mat');
data = S.datasa;

X = flip(data.v1,2)/1000;    % [mm]
Z = flip(data.v3,2)/1000;    % [mm]
u = flip(data.v4,2);
w = -flip(data.v6,2);
dx = abs(X(2,1) - X(1,1));   
dz = abs(Z(1,2) - Z(1,1)); 
[du_dz, du_dx] = gradient(u, dx, dz);
[dw_dz, dw_dx] = gradient(w, dx, dz);

curl = dw_dx - du_dz;               
div_manual = du_dx + dw_dz;

div2D = data.v21;

fprintf('max abs difference = %g\n', max(abs(div2D(:) - div_manual(:)),[],'omitnan'));