function [p, Mask] = directIntegrateNS2D(X, Z, u2, w2, ut, wt, ...
    rho, mu, valid, seedMask, pSeedFull, weightingMode)
%DIRECTINTEGRATENS2D Strong-form marching integration using unsteady NS.
%
% Inputs:
%   X,Z        [m] same size as velocities
%   u1,w1      at t1
%   u2,w2      at t2 (target time)
%   u3,w3      at t3
%   rho, mu    density [kg/m^3], viscosity [Pa*s] (set mu=0 for inviscid)
%   valid      logical mask of usable PIV nodes
%   seedMask   logical mask of nodes where pressure is known (ROI seeds)
%   pSeedFull  full-size array with pressure on seedMask, NaN elsewhere
%   weightingMode: 'equal' or 'isopotential' (uses U,W direction like your old code)
%
% Outputs:
%   p          pressure on valid nodes (NaN outside valid)
%   Mask       -1 outside valid, 0 unknown, 1 known

if nargin < 17 || isempty(weightingMode)
    weightingMode = 'equal';
end
weightingMode = lower(weightingMode);

% --- init Mask like your old flow.Mask convention
Mask = -1 * ones(size(X));
Mask(valid) = 0;

% seeds = known
Mask(seedMask & valid) = 1;

% init pressure
p = nan(size(X));
p(seedMask & valid) = pSeedFull(seedMask & valid);

p_inf = 102197;   % Pa
V_inf = 15;       % m/s
q_inf = 0.5 * rho * V_inf^2;
Vsq = u2.^2 + w2.^2;

% --- compute grid steps (assume rectilinear uniform)
dx = median(abs(diff(X(:,1))), 'omitnan');
dz = median(abs(diff(Z(1,:))), 'omitnan');


% spatial derivatives at t2
[uz, ux] = gradient(u2, dz, dx);
[wz, wx] = gradient(w2, dz, dx);

% convective accel
ax = ut + u2.*ux + w2.*uz;
az = wt + u2.*wx + w2.*wz;

% viscous term (optional)
if mu ~= 0
    [uxx, ~ ] = gradient(ux, dx, dz);
    [~ , uzz] = gradient(uz, dx, dz);
    [wxx, ~ ] = gradient(wx, dx, dz);
    [~ , wzz] = gradient(wz, dx, dz);
    lap_u = uxx + uzz;
    lap_w = wxx + wzz;
else
    lap_u = 0;
    lap_w = 0;
end

% pressure gradients from momentum
PX = -rho*ax + mu*lap_u;   % dp/dx
PZ = -rho*az + mu*lap_w;   % dp/dz

% kill gradients outside valid
PX(~valid) = NaN;
PZ(~valid) = NaN;

% --- build initial boundary: unknown nodes adjacent to known nodes
kernel = ones(3,3); kernel(2,2) = 0;
isBoundary = (Mask==0) & (convn(Mask==1, kernel, 'same') > 0);
minP = min(min(p));
maxP = max(max(p));
Boundary = find(isBoundary);
itNum=0;
% --- marching loop
while any(Mask(:)==0)
    itNum=itNum+1;
    if isempty(Boundary)
        warning('Marching stalled: no boundary candidates but unknown nodes remain.');
        break;
    end

    for idxLin = Boundary(:).'
        [i,k] = ind2sub(size(Mask), idxLin);

        % skip if became known earlier this sweep
        if Mask(i,k) ~= 0, continue; end

        % compute p(i,k) from known neighbors
        p_here = calcP_fromNeighbors_NS2D(X,Z,u2,w2, p, Mask, PX,PZ, i,k, weightingMode);

        if isfinite(p_here)
            p(i,k) = p_here;
            Mask(i,k) = 2; % temporary known
        end
    end

    if (rem(itNum,20)==0)
        % plotNormal(X,Z,Mask,"Mask",1,-1.5,1,1.5)
        % plotNormal(X,Z,grid.insideMask,"insideMask",2,-1.5,1,1.5)
        cp    = (p - p_inf) / q_inf;
        cpbern = 1 - Vsq / V_inf^2;

        plotCompare(X,Z,cpbern,cp,"Bernoulli cp", "Integrated cp",10,-0.8,0.1,1)
        %plotNormal(X,Z,cp, "Integrated Pressure",10,-0.8,0.1,1)
    end
    Mask(Mask==2) = 1;

    % new boundary
    isBoundary = (Mask==0) & (convn(Mask==1, kernel, 'same') > 0);
    Boundary = find(isBoundary);
end

% outside valid -> NaN
p(~valid) = NaN;

end