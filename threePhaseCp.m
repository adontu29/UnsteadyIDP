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

% Rectangle ROI in mm: [xL xR zB zT]
rect_mm = [-80, 5, -87, 87];
rect_m  = rect_mm / 1000;

rho = 1.225;  % kg/m^3

theta = [0.83, 12.435, 31.0340]*pi/180;   % parameter values for the 3 snapshots
fHz   = [1.0, 1.0, 1.0]*2.39;             % (your choice)

% If theta is PHASE (rad), then dtheta/dt = 2*pi*f
dtheta_dt = 2*pi*fHz;     % rad/s

% ------------------------------------------------------------
% Preallocate storage for full fields (needed for du/dt, dw/dt)
% ------------------------------------------------------------
outAll   = cell(3,1);
phiROI   = cell(3,1);

Xall     = cell(3,1);   % [m]
Zall     = cell(3,1);   % [m]
Uall     = cell(3,1);   % [m/s]
Wall     = cell(3,1);   % [m/s]
validAll = cell(3,1);   % logical

% ------------------------------------------------------------
% Loop over snapshots: build masks, store full fields, solve phi in ROI
% ------------------------------------------------------------
for k = 1:3
    S = load(fullfile(dataDir, fileList(k)));
    d = S.datasa;

    % Flip as you do
    Xmm = flip(d.v1,2);    % mm
    Zmm = flip(d.v3,2);    % mm
    u   = flip(d.v4,2);
    w   = -flip(d.v6,2);

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

    % Clean NaNs in velocity (keep mask as authority)
    u(~isfinite(u)) = 0;
    w(~isfinite(w)) = 0;

    if k==2
        figure;
        plotNormal(Xmm, Zmm, valid, "Mask", 1, -0.5, 1, 1.5);
        title("Valid mask (snapshot 2)");
    end

    % Store FULL fields (meters for X/Z)
    Xall{k}     = Xmm/1000;   % m
    Zall{k}     = Zmm/1000;   % m
    Uall{k}     = u;
    Wall{k}     = w;
    validAll{k} = valid;

    % Solve Laplace potential in ROI (inputs in meters)
    out = solveLaplacePotentialRect(Xmm/1000, Zmm/1000, -u, -w, valid, rect_m);

    outAll{k} = out;
    phiROI{k} = -out.phi_roi;   % your sign convention
end

% Remove arbitrary constant in phi inside ROI for each snapshot
for k=1:3
    tmp = phiROI{k};
    tmp = tmp - mean(tmp(outAll{k}.mask_roi), 'all', 'omitnan');
    phiROI{k} = tmp;
end

% ------------------------------------------------------------
% Unequal-theta 3-point weights at snapshot 2
% ------------------------------------------------------------
th1 = theta(1); th2 = theta(2); th3 = theta(3);

weight1 = (th2 - th3)/((th1 - th2)*(th1 - th3));
weight2 = (2*th2 - th1 - th3)/((th2 - th1)*(th2 - th3));
weight3 = (th2 - th1)/((th3 - th1)*(th3 - th2));

% ------------------------------------------------------------
% dphi/dt at snapshot 2 (ROI)
% ------------------------------------------------------------
phi1 = phiROI{1}; phi2 = phiROI{2}; phi3 = phiROI{3};

dphi_dtheta_2 = weight1*phi1 + weight2*phi2 + weight3*phi3;
dphi_dt_2     = dphi_dtheta_2 * dtheta_dt(2);

% ------------------------------------------------------------
% Pull full fields at snapshot 2 (and 1 & 3 for time derivatives)
% ------------------------------------------------------------
out2   = outAll{2};
i1 = out2.idx(1); i2 = out2.idx(2);
j1 = out2.idx(3); j2 = out2.idx(4);

X2 = Xall{2};  Z2 = Zall{2};      % m
u1 = Uall{1};  w1 = Wall{1};
u2 = Uall{2};  w2 = Wall{2};
u3 = Uall{3};  w3 = Wall{3};
valid2 = validAll{2};

% ------------------------------------------------------------
% Time derivatives ut, wt at snapshot 2 using unequal-theta stencil
% ------------------------------------------------------------
du_dtheta_2 = weight1*u1 + weight2*u2 + weight3*u3;
dw_dtheta_2 = weight1*w1 + weight2*w2 + weight3*w3;

ut_2 = du_dtheta_2 * dtheta_dt(2);
wt_2 = dw_dtheta_2 * dtheta_dt(2);

% ------------------------------------------------------------
% Bernoulli pressure in ROI at snapshot 2
% ------------------------------------------------------------
% Velocity in ROI
uROI = u2(i1:i2, j1:j2);
wROI = w2(i1:i2, j1:j2);
V2sq = uROI.^2 + wROI.^2;

maskROI = out2.mask_roi;

% Potential velocity from phi (ensure out2.dx/out2.dz are meters)
[dphi_di, dphi_dj] = gradient(out2.phi_roi, out2.dx, out2.dz);

% NOTE: these signs depend on your phi convention; keep as you used
u_phi = uROI;
w_phi = wROI;
Vphisq = u_phi.^2 + w_phi.^2;

qTerm = -rho*dphi_dt_2 - 0.5*rho*Vphisq;

p_inf = 102197;   % Pa
V_inf = 15;       % m/s
q_inf = 0.5 * rho * V_inf^2;

% Reference using left edge of ROI
nc = size(qTerm,2);
Nleft = max(2, round(0.05*nc));

leftMask = false(size(qTerm));
leftMask(:,1:Nleft) = true;
leftMask = leftMask & maskROI & isfinite(qTerm);

qLeftMean = mean(qTerm(leftMask), 'all', 'omitnan');

pROI_abs = qTerm - qLeftMean + p_inf;      % Pa
CpROI    = (pROI_abs - p_inf) / q_inf;

% ------------------------------------------------------------
% Quick plots (ROI)
% ------------------------------------------------------------
CpPlot = CpROI;
CpPlot(~maskROI) = NaN;

minCp = min(CpPlot, [], 'all', 'omitnan');
maxCp = max(CpPlot, [], 'all', 'omitnan');

figure(1)
hold on
plotNormal(X2(i1:i2,j1:j2)*1000, Z2(i1:i2,j1:j2)*1000, CpPlot, ...
    "Cp", 1, minCp, (maxCp-minCp)/10, maxCp);
title('Cp in upstream ROI (snapshot 2)');
xlabel('z [mm]'); ylabel('x [mm]');

Cp_speed = 1 - V2sq / V_inf^2;
Cp_speed(~maskROI) = NaN;
figure;
plotNormal(X2(i1:i2,j1:j2)*1000, Z2(i1:i2,j1:j2)*1000, Cp_speed, ...
    "Cp\_speed", 2, -1, 0.1, 1);
title('Cp from speed only (sanity check)');

% ------------------------------------------------------------
% Build full-grid seed mask + seed pressure field (NOW that pROI_abs exists)
% ------------------------------------------------------------
seedMaskFull = false(size(u2));
seedMaskFull(i1:i2, j1:j2) = maskROI;

pSeedFull = nan(size(u2));
pSeedFull(i1:i2, j1:j2) = pROI_abs;

% ------------------------------------------------------------
% Marching integration using NS gradients
% ------------------------------------------------------------
% weightingMode can be "equal" or "isopotential" depending on your marcher
weightingMode = "equal";

% NOTE: your directIntegrateNS2D must accept these inputs:
% (X,Z,u,w,ut,wt,rho,mu,valid,seedMask,pSeed,weightingMode)
[p, IntMask] = directIntegrateNS2D( ...
    X2, Z2, u2, w2, ut_2, wt_2, ...
    rho, 0, valid2, seedMaskFull, pSeedFull, weightingMode);
