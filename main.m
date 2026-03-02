clear
close all

S    = load('PIVdata/U15A10-10K100P0000.mat');
data = S.datasa;

X = data.v1;
Z = data.v3;
u = data.v4;
w = data.v6;

dx = abs(X(2,1) - X(1,1))/1000;    
dz = abs(Z(1,2) - Z(1,1))/1000;  

[du_dx, ~] = gradient(u, dx, dz);
[~, dw_dz] = gradient(w, dx, dz);

% Poisson RHS (or use zeros for Laplace)


omegaMag = data.v20;
div2D    = data.v21;
isValid  = data.v44;

f =div2D;
 %f = du_dx + dw_dz;
% f = zeros(size(div2D));   % if you really want Laplace

valid    = (isValid == 1);

% Just for *using* the solution later (not for domain)
omega_thr = prctile(abs(omegaMag(valid)), 40);
div_thr   = prctile(abs(div2D(valid)),   80);
potMask   = valid & abs(omegaMag) <= omega_thr & abs(div2D) <= div_thr;


% --- keep only the largest connected component ---
cc = bwconncomp(valid, 4);             % 4-connectivity
numPixels = cellfun(@numel, cc.PixelIdxList);
[~, iBig] = max(numPixels);

validBig = false(size(valid));
validBig(cc.PixelIdxList{iBig}) = true;
valid = validBig;


% --- erode the outer boundary by n cells ---
nErode = 3;   % 2 or 3, tweak as you like

% Option A: bwmorph (iterated erosion)
valid = bwmorph(valid, 'erode', nErode);


fprintf('NaNs in u over valid: %d\n', nnz(valid & isnan(u)));
fprintf('NaNs in w over valid: %d\n', nnz(valid & isnan(w)));
fprintf('NaNs in f over valid: %d\n', nnz(valid & isnan(f)));

f(~valid) = 0;
u(~isfinite(u)) = 0;
w(~isfinite(w)) = 0;
f(~isfinite(f)) = 0;


phi = solvePoissonIrregularNeumann(X, Z, f, valid, u, w, 0);
% Then call the solver with validClean instead of validMask


% Visualise: only trust phi in potMask
figure;
contourf(X, Z, phi, 30); axis equal tight; colorbar;
figure
contour(X,Z,valid,[1 1],'k','LineWidth',1);    % measurement domain
figure
contour(X,Z,potMask,[1 1],'r','LineWidth',1.5);% potential region
title('\phi, valid domain (black), potential region (red)');

% Potential velocity
[phi_x, phi_z] = gradient(phi, dx, dz);  

[phi_x, phi_z] = gradient(phi, dx, dz);

fprintf('max |u(valid)|     = %e\n', max(abs(u(valid)), [], 'all'));
fprintf('max |phi_x(valid)| = %e\n', max(abs(phi_x(valid)), [], 'all'));
fprintf('max |w(valid)|     = %e\n', max(abs(w(valid)), [], 'all'));
fprintf('max |phi_z(valid)| = %e\n', max(abs(phi_z(valid)), [], 'all'));

mask = potMask & ~isnan(u) & ~isnan(w);
up = phi_x;   wp = phi_z;

e_u  = u(mask) - up(mask);
e_w  = w(mask) - wp(mask);

fprintf('RMS error u: %g m/s\n', rms(e_u));
fprintf('RMS error w: %g m/s\n', rms(e_w));