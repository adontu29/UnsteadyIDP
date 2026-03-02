S    = load('U15A10-10K100P0000.mat');
data = S.datasa;

X = data.v1;
Z = data.v3;
u = data.v4;
w = data.v6;
dx = abs(X(2,1) - X(1,1))/1000;   
dz = abs(Z(1,2) - Z(1,1))/1000; 
[du_dx, ~] = gradient(u, dx, dz);
[~, dw_dz] = gradient(w, dx, dz);
div_manual = du_dx + dw_dz;

div2D = data.v21;

fprintf('max abs difference = %g\n', max(abs(div2D(:) - div_manual(:)),[],'omitnan'));