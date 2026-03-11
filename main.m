clear
close all

S    = load('PIVdata/U15A10-10K100P0000.mat');
data = S.datasa;


X = flip(data.v1,2);    % [mm]
Z = flip(data.v3,2);    % [mm]
u = flip(data.v4,2);
w = -flip(data.v6,2);

dx = abs(X(2,1) - X(1,1))/1000;    
dz = abs(Z(1,2) - Z(1,1))/1000;  

[du_dz, ~] = gradient(u, dz, dx);
[~, dw_dx] = gradient(w, dz, dx);

% Poisson RHS (or use zeros for Laplace)


omegaMag = data.v20;
div2D    = flip(data.v21,2);
isValid  = flip(data.v44,2);

f =div2D;
%f = du_dx + dw_dz;
% f = zeros(size(div2D));   % if you really want Laplace

valid    = (isValid == 1);

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

rect = [-80, 5,-87, 87 ];


plotNormal(X,Z,valid,"Mask",1,-0.5,1,1.5)

xL = -102;
xR = -20;
zB = -87;
zT = 87;

out = solveLaplacePotentialRect(X/1000, Z/1000, u, w, valid, rect/1000);
phi = out.phi;

minPhi = min(min(phi));
maxPhi = max(max(phi));

figure(1)
hold on
plotNormal(X,Z,phi,"phi",1,minPhi,(maxPhi-minPhi)/10,maxPhi)



[phi_z, phi_x] = gradient(-out.phi_roi, out.dz, out.dx);


alpha_phi = atan2(mean(phi_x(out.mask_roi),'all'), mean(phi_z(out.mask_roi),'all')) * 180/pi

err_u = u(out.idx(1):out.idx(2), out.idx(3):out.idx(4)) - phi_x;
err_w = w(out.idx(1):out.idx(2), out.idx(3):out.idx(4)) - phi_z;

figure
quiver(X(out.idx(1):out.idx(2), out.idx(3):out.idx(4)),Z(out.idx(1):out.idx(2), out.idx(3):out.idx(4)),u(out.idx(1):out.idx(2), out.idx(3):out.idx(4)),w(out.idx(1):out.idx(2), out.idx(3):out.idx(4)), 0.5)
title("u,v")

figure
quiver(X(out.idx(1):out.idx(2), out.idx(3):out.idx(4)),Z(out.idx(1):out.idx(2), out.idx(3):out.idx(4)),phi_x,phi_z, 0.5)
title("phi_x,phi_z")

figure
quiver(X(out.idx(1):out.idx(2), out.idx(3):out.idx(4)),Z(out.idx(1):out.idx(2), out.idx(3):out.idx(4)),err_u,err_w, 1)
title("error plot")


fprintf('RMS u error: %.3g\n', rms(err_u(out.mask_roi), 'omitnan'));
fprintf('RMS w error: %.3g\n', rms(err_w(out.mask_roi), 'omitnan'));