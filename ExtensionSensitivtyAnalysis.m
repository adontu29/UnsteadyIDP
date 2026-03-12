clear; close all; clc;

% ------------------------------------------------------------
% Load one snapshot for the sensitivity study
% ------------------------------------------------------------
dataDir = "PIVdata";
fileName = "U15A10-10K100P0150.mat";   % snapshot 2 for example

S = load(fullfile(dataDir, fileName));
d = S.datasa;

X = flip(d.v1,2) / 1000;   % m
Z = flip(d.v3,2) / 1000;   % m
u = flip(d.v4,2);
w = -flip(d.v6,2);

isValid = flip(d.v44,2);
valid = (isValid == 1);

% Largest connected component
cc = bwconncomp(valid, 4);
numPixels = cellfun(@numel, cc.PixelIdxList);
[~, iBig] = max(numPixels);
validBig = false(size(valid));
validBig(cc.PixelIdxList{iBig}) = true;
valid = validBig;

% Erode mask without shrinking holes
nErode = 3;
holes = imfill(valid, 'holes') & ~valid;
validFilled = imfill(valid, 'holes');
validErodedFilled = bwmorph(validFilled, 'erode', nErode);
valid = validErodedFilled & ~holes;

u(~isfinite(u)) = 0;
w(~isfinite(w)) = 0;

% ------------------------------------------------------------
% Kept ROI
% ------------------------------------------------------------
rectKeep_mm = [-80, 5, -87, 87];
rectKeep_m  = rectKeep_mm / 1000;

% ------------------------------------------------------------
% Freestream in lab frame and pitch angle
% ------------------------------------------------------------
U_inf_lab = 15;
theta_deg = 12.435;
theta = theta_deg * pi/180;

% Rotate freestream into airfoil coordinates
Uinf_body =  U_inf_lab * cos(theta);
Winf_body = -U_inf_lab * sin(theta);   % flip sign if your convention needs it

% ------------------------------------------------------------
% Extension lengths to test [m]
% These are EXTRA lengths added upstream beyond xL_keep
% ------------------------------------------------------------
extLengths_mm = [0 5 10 20 30 40 60 80];
extLengths_m  = extLengths_mm / 1000;

% Left boundary of kept ROI
xL_keep = rectKeep_m(1);

% ------------------------------------------------------------
% Storage
% ------------------------------------------------------------
nCases = numel(extLengths_m);

rmse_u   = nan(nCases,1);
rmse_w   = nan(nCases,1);
rmse_vec = nan(nCases,1);
mean_ang = nan(nCases,1);

outCases = cell(nCases,1);

% ------------------------------------------------------------
% Loop over extension lengths
% ------------------------------------------------------------
for k = 1:nCases

    Lext = extLengths_m(k);

    rectExt_m = rectKeep_m;
    rectExt_m(1) = xL_keep - Lext;   % move left boundary upstream

    out = solveLaplacePotentialRectExtended( ...
        X, Z, u, w, valid, ...
        rectExt_m, rectKeep_m, ...
        Uinf_body, Winf_body);

    outCases{k} = out;

    % --------------------------------------------------------
    % Compute reconstructed velocity from phi in kept ROI
    %
    % IMPORTANT:
    % You said this swap is intentional for your PIV indexing.
    % Keep it exactly like this.
    % --------------------------------------------------------
    [phi_z, phi_x] = gradient(-out.phi_keep, out.dz, out.dx);

    % Measured velocities in kept ROI
    i1 = out.keep_idx(1);
    i2 = out.keep_idx(2);
    j1 = out.keep_idx(3);
    j2 = out.keep_idx(4);

    uROI = u(i1:i2, j1:j2);
    wROI = w(i1:i2, j1:j2);

    mask = out.mask_keep;

    % Errors
    err_u = uROI - phi_x;
    err_w = wROI - phi_z;

    % Apply mask
    err_u(~mask) = NaN;
    err_w(~mask) = NaN;

    % RMSEs
    rmse_u(k) = sqrt(mean(err_u(mask).^2, 'omitnan'));
    rmse_w(k) = sqrt(mean(err_w(mask).^2, 'omitnan'));
    rmse_vec(k) = sqrt(mean(err_u(mask).^2 + err_w(mask).^2, 'omitnan'));

    % Mean reconstructed flow angle (optional diagnostic)
    mean_ang(k) = atan2(mean(phi_x(mask), 'omitnan'), mean(phi_z(mask), 'omitnan')) * 180/pi;

    fprintf('Extension = %6.1f mm | RMSE_u = %.4f | RMSE_w = %.4f | RMSE_vec = %.4f\n', ...
        extLengths_mm(k), rmse_u(k), rmse_w(k), rmse_vec(k));
end

% ------------------------------------------------------------
% Plot sensitivity
% ------------------------------------------------------------
figure;
plot(extLengths_mm, rmse_u, '-o', 'LineWidth', 1.5);
hold on
plot(extLengths_mm, rmse_w, '-s', 'LineWidth', 1.5);
plot(extLengths_mm, rmse_vec, '-d', 'LineWidth', 1.5);
grid on
xlabel('Upstream extension length [mm]');
ylabel('RMSE [m/s]');
legend('RMSE_u', 'RMSE_w', 'RMSE_{vec}', 'Location', 'best');
title('Sensitivity of velocity reconstruction to extension length');

figure;
plot(extLengths_mm, mean_ang, '-o', 'LineWidth', 1.5);
grid on
xlabel('Upstream extension length [mm]');
ylabel('Mean reconstructed angle [deg]');
title('Mean angle from reconstructed potential velocity');