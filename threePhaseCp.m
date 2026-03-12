clear; close all; clc;

% ------------------------------------------------------------
% User inputs
% ------------------------------------------------------------
dataDir = "PIVdata";
fileList = [
    "U15A10-10K100P0000.mat"
    "U15A10-10K100P0150.mat"
    "U15A10-10K100P0300.mat"
];

% Original ROI that you want to KEEP in the end
rectKeep_mm = [-80, 5, -87, 87];
rectKeep_m  = rectKeep_mm / 1000;

% Extended ROI used only for the Laplace solve
% Change only the LEFT boundary
rectExt_mm = [-80, 5, -87, 87];
rectExt_m  = rectExt_mm / 1000;

rho = 1.225;      % kg/m^3
p_inf = 102197;   % Pa
U_inf_lab = 15;   % m/s (known in lab frame)

theta = [0.83, 12.435, 31.0340] * pi/180;   % pitch angle for each snapshot
fHz   = [1.0, 1.0, 1.0] * 2.39;
dtheta_dt = 2*pi*fHz;   % rad/s

% ------------------------------------------------------------
% Storage
% ------------------------------------------------------------
outAll    = cell(3,1);
phiExtAll = cell(3,1);   % phi on extended domain
phiKeepAll = cell(3,1);  % phi only on original keep-ROI

Xall     = cell(3,1);    % full measured grid [m]
Zall     = cell(3,1);
Uall     = cell(3,1);
Wall     = cell(3,1);
validAll = cell(3,1);

% ------------------------------------------------------------
% Loop over snapshots
% ------------------------------------------------------------
for k = 1:3
    S = load(fullfile(dataDir, fileList(k)));
    d = S.datasa;

    Xmm = flip(d.v1,2);       % mm
    Zmm = flip(d.v3,2);       % mm
    u   = flip(d.v4,2);       % airfoil-frame x velocity
    w   = -flip(d.v6,2);      % airfoil-frame z velocity (your current convention)

    isValid = flip(d.v44,2);
    valid = (isValid == 1);

    % Keep largest connected component
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

    % Clean NaNs in velocity (mask remains the authority)
    u(~isfinite(u)) = 0;
    w(~isfinite(w)) = 0;

    if k == 2
        figure;
        plotNormal(Xmm, Zmm, valid, "Mask", 1, -0.5, 1, 1.5);
        title("Valid mask (snapshot 2)");
    end

    % Store full measured fields
    Xall{k}     = Xmm / 1000;
    Zall{k}     = Zmm / 1000;
    Uall{k}     = u;
    Wall{k}     = w;
    validAll{k} = valid;

    % --------------------------------------------------------
    % Rotate freestream from lab frame to AIRFOIL coordinates
    %
    % Lab freestream is [U_inf_lab; 0].
    % In airfoil/body coordinates:
    %   Ux = U_inf_lab*cos(theta)
    %   Uz = -U_inf_lab*sin(theta)
    %
    % IMPORTANT:
    % If your z-axis/sign convention is opposite, flip the sign of Winf_body.
    % --------------------------------------------------------
    Uinf_body =  U_inf_lab * cos(theta(k));
    Winf_body =  U_inf_lab * sin(theta(k));

    % --------------------------------------------------------
    % Solve Laplace on extended domain
    % --------------------------------------------------------
    out = solveLaplacePotentialRectExtended( ...
        Xmm/1000, Zmm/1000, ...
        u, w, valid, ...
        rectExt_m, rectKeep_m, ...
        Uinf_body, Winf_body);

    outAll{k}     = out;
    phiExtAll{k}  = out.phi_ext;
    phiKeepAll{k} = out.phi_keep;

    if k == 2
        figure;
        phiPlot = out.phi_ext;
        phiPlot(~out.mask_ext) = NaN;
        imagesc(out.Zext(1,:)*1000, out.Xext(:,1)*1000, phiPlot);
        axis xy equal tight
        colorbar
        xlabel('z [mm]');
        ylabel('x [mm]');
        title('Extended-domain potential \phi (snapshot 2)');
    end
end

% ------------------------------------------------------------
% Unequal-theta 3-point weights at snapshot 2
% ------------------------------------------------------------
th1 = theta(1); th2 = theta(2); th3 = theta(3);

weight1 = (th2 - th3)/((th1 - th2)*(th1 - th3));
weight2 = (2*th2 - th1 - th3)/((th2 - th1)*(th2 - th3));
weight3 = (th2 - th1)/((th3 - th1)*(th3 - th2));

% ------------------------------------------------------------
% dphi'/dt on the EXTENDED domain
% ------------------------------------------------------------
phiExt1 = phiExtAll{1};
phiExt2 = phiExtAll{2};
phiExt3 = phiExtAll{3};

dphi_dtheta_ext_2 = weight1*phiExt1 + weight2*phiExt2 + weight3*phiExt3;
dphi_dt_ext_2     = dphi_dtheta_ext_2 * dtheta_dt(2);

% ------------------------------------------------------------
% Estimate c_dot from the artificial upstream strip
% ------------------------------------------------------------
out2 = outAll{2};
refMaskExt = out2.ref_mask_ext;

c_dot_2 = mean(dphi_dt_ext_2(refMaskExt), 'all', 'omitnan');

% Corrected time derivative
dphi_dt_ext_corr_2 = dphi_dt_ext_2 - c_dot_2;

% Keep only original measured ROI part
dphi_dt_keep_corr_2 = dphi_dt_ext_corr_2(out2.keep_rows, :);

% Also keep phi itself in the original ROI, if needed later
phi2_keep = out2.phi_keep;
maskKeep  = out2.mask_keep;

% ------------------------------------------------------------
% Pull measured fields at snapshot 2
% ------------------------------------------------------------
X2 = Xall{2};
Z2 = Zall{2};

u1 = Uall{1};  w1 = Wall{1};
u2 = Uall{2};  w2 = Wall{2};
u3 = Uall{3};  w3 = Wall{3};

valid2 = validAll{2};

i1 = out2.keep_idx(1);
i2 = out2.keep_idx(2);
j1 = out2.keep_idx(3);
j2 = out2.keep_idx(4);

% ------------------------------------------------------------
% Time derivatives ut, wt at snapshot 2 on FULL measured grid
% ------------------------------------------------------------
du_dtheta_2 = weight1*u1 + weight2*u2 + weight3*u3;
dw_dtheta_2 = weight1*w1 + weight2*w2 + weight3*w3;

ut_2 = du_dtheta_2 * dtheta_dt(2);
wt_2 = dw_dtheta_2 * dtheta_dt(2);

% ------------------------------------------------------------
% Bernoulli pressure in kept ROI
%
% p = p_inf + 0.5*rho*(U_inf^2 - V^2) - rho*dphi/dt_corrected
%
% Here U_inf is the MAGNITUDE of the lab freestream (= 15 m/s),
% which is invariant under rotation.
% ------------------------------------------------------------
uROI = u2(i1:i2, j1:j2);
wROI = w2(i1:i2, j1:j2);
V2sq = uROI.^2 + wROI.^2;

pROI_abs = p_inf + 0.5*rho*(U_inf_lab^2 - V2sq) - rho*dphi_dt_keep_corr_2;
pROI_abs(~maskKeep) = NaN;

q_inf = 0.5 * rho * U_inf_lab^2;
CpROI = (pROI_abs - p_inf) / q_inf;
CpROI(~maskKeep) = NaN;

% ------------------------------------------------------------
% Quick plots
% ------------------------------------------------------------
minCp = min(CpROI, [], 'all', 'omitnan');
maxCp = max(CpROI, [], 'all', 'omitnan');

figure;
plotNormal(X2(i1:i2,j1:j2)*1000, Z2(i1:i2,j1:j2)*1000, CpROI, ...
    "Cp", 1, minCp, (maxCp-minCp)/10, maxCp);
title('Cp in kept ROI (snapshot 2)');
xlabel('z [mm]');
ylabel('x [mm]');

Cp_speed = 1 - V2sq / U_inf_lab^2;
Cp_speed(~maskKeep) = NaN;

figure;
plotNormal(X2(i1:i2,j1:j2)*1000, Z2(i1:i2,j1:j2)*1000, Cp_speed, ...
    "Cp\_speed", 2, -1, 0.1, 1);
title('Cp from speed only (sanity check)');

% ------------------------------------------------------------
% Build full-grid seed mask + seed pressure field
% ------------------------------------------------------------
seedMaskFull = false(size(u2));
seedMaskFull(i1:i2, j1:j2) = maskKeep;

pSeedFull = nan(size(u2));
pSeedFull(i1:i2, j1:j2) = pROI_abs;

% ------------------------------------------------------------
% Marching integration using NS gradients
% ------------------------------------------------------------
weightingMode = "equal";

[p, IntMask] = directIntegrateNS2D( ...
    X2, Z2, u2, w2, ut_2, wt_2, ...
    rho, 0, valid2, seedMaskFull, pSeedFull, weightingMode);