function out = solveLaplacePotentialRect(X, Z, u, w, valid, rect)
%SOLVELAPLACEPOTENTIALRECT Solve Laplace for potential phi in a rectangular ROI
%   ∇²phi = 0 on ROI
%   Neumann BC from PIV: dphi/dx = u, dphi/dz = w on boundary
%
% Inputs
%   X,Z    : grids (same size) [mm or any units, but consistent]
%   u,w    : velocities (same size) [units per second]; interpreted as:
%            u = dphi/dx, w = dphi/dz  
%   valid  : logical mask (same size) selecting usable data
%   rect   : [x_left x_right z_bottom z_top] in same units as X,Z
%
% Output struct 'out' fields:
%   out.phi       : phi on full grid (NaN outside ROI/mask)
%   out.phi_roi   : phi on ROI grid only
%   out.mask_roi  : mask used inside ROI
%   out.idx       : [i1 i2 j1 j2] index bounds (rows, cols)
%   out.rect      : rect used
%
% Notes
%   - Pure Neumann problem has an arbitrary constant; we fix phi(ref)=0.
%   - Assumes Cartesian grid with uniform spacing in x and z.


    arguments
        X (:,:) double
        Z (:,:) double
        u (:,:) double
        w (:,:) double
        valid (:,:) logical
        rect (1,4) double
    end

    % --- unpack rectangle bounds (allow user to give swapped order) ---
    xL = min(rect(1), rect(2));
    xR = max(rect(1), rect(2));
    zB = min(rect(3), rect(4));
    zT = max(rect(3), rect(4));

    % --- find index bounds for ROI using X as "x", Z as "z" ---
    % X varies mainly along rows (i), Z varies mainly along cols (j) in your data.
    xline = X(:,1);
    zline = Z(1,:);

    i1 = find(xline >= xL, 1, "first");
    i2 = find(xline <= xR, 1, 'last');
    j1 = find(zline >= zB, 1, 'first');
    j2 = find(zline <= zT, 1, 'last');

    if isempty(i1) || isempty(i2) || isempty(j1) || isempty(j2) || i2<i1 || j2<j1
        error('ROI bounds not found on grid. Check rect and X/Z axes directions.');
    end

    % ROI slices
    Xr = X(i1:i2, j1:j2);
    Zr = Z(i1:i2, j1:j2);
    ur = u(i1:i2, j1:j2);
    wr = w(i1:i2, j1:j2);
    mr = valid(i1:i2, j1:j2);

    % --- estimate uniform spacings in METERS (or in same units consistently) ---

    dx = median(abs(diff(Xr(:,1))));
    dz = median(abs(diff(Zr(1,:))));

    if ~isfinite(dx) || ~isfinite(dz) || dx<=0 || dz<=0
        error('Bad dx/dz computed. Check X/Z in ROI.');
    end

    % --- build unknown list only where mask is true ---
    [nr, nc] = size(Xr);
    id = zeros(nr, nc);              % map (i,j) -> unknown index
    lin = find(mr);
    N = numel(lin);
    if N < 10
        error('Too few valid points inside ROI to solve (%d).', N);
    end
    id(lin) = 1:N;

    % sparse system
    A = spalloc(N, N, 5*N);
    b = zeros(N,1);

    idx = @(ii,jj) id(ii,jj);

    % helper: safe access to boundary velocity (use local node value)
    % u = dphi/dx, w = dphi/dz
    for ii = 1:nr
        for jj = 1:nc
            p = idx(ii,jj);
            if p==0, continue; end

            % Base Laplace stencil coefficients for this node:
            % (phiE - 2phiP + phiW)/dx^2 + (phiN - 2phiP + phiS)/dz^2 = 0
            cP = -2*(1/dx^2 + 1/dz^2);

            % --- West neighbor (ii-1, jj) ---
            if ii-1 >= 1 && idx(ii-1,jj) ~= 0
                q = idx(ii-1,jj);
                A(p,q) = A(p,q) + 1/dx^2;
            else
                % Neumann on west boundary: u = dphi/dx
                % ghost phiW = phiP - dx*u  (first-order)
                % contributes: (phiE - phiP - dx*u)/dx^2  => RHS += -u/dx
                b(p) = b(p) - ur(ii,jj)/dx;
                % no neighbor coefficient added, and P coefficient changes from -2/dx^2 to -1/dx^2
                cP = cP + 1/dx^2;
            end

            % --- East neighbor (ii+1, jj) ---
            if ii+1 <= nr && idx(ii+1,jj) ~= 0
                q = idx(ii+1,jj);
                A(p,q) = A(p,q) + 1/dx^2;
            else
                % Neumann east: ghost phiE = phiP + dx*u
                % contributes: (phiW - phiP + dx*u)/dx^2 => RHS += +u/dx
                b(p) = b(p) + ur(ii,jj)/dx;
                cP = cP + 1/dx^2;
            end

            % --- South neighbor (ii, jj-1)  (z decreases downward if your zline increases rightward)
            if jj-1 >= 1 && idx(ii,jj-1) ~= 0
                q = idx(ii,jj-1);
                A(p,q) = A(p,q) + 1/dz^2;
            else
                % Neumann south (bottom): ghost phiS = phiP - dz*w
                % term becomes (phiN - phiP - dz*w)/dz^2 => RHS += -w/dz
                b(p) = b(p) - wr(ii,jj)/dz;
                cP = cP + 1/dz^2;
            end

            % --- North neighbor (ii, jj+1) ---
            if jj+1 <= nc && idx(ii,jj+1) ~= 0
                q = idx(ii,jj+1);
                A(p,q) = A(p,q) + 1/dz^2;
            else
                % Neumann north (top): ghost phiN = phiP + dz*w
                % term becomes (phiS - phiP + dz*w)/dz^2 => RHS += +w/dz
                b(p) = b(p) + wr(ii,jj)/dz;
                cP = cP + 1/dz^2;
            end

            A(p,p) = A(p,p) + cP;
        end
    end

    % --- Fix the nullspace of pure Neumann: set one reference point phi=0 ---
    ref = 1;              % first unknown
    A(ref,:) = 0;
    A(ref,ref) = 1;
    b(ref) = 0;

    % --- solve ---
    phi_vec = A \ b;

    % --- pack back into ROI grid and full grid ---
    phi_roi = nan(nr,nc);
    phi_roi(lin) = phi_vec;

    phi_full = nan(size(X));
    phi_full(i1:i2, j1:j2) = phi_roi;

    out = struct();
    out.phi      = phi_full;
    out.phi_roi  = phi_roi;
    out.mask_roi = mr;
    out.idx      = [i1 i2 j1 j2];
    out.rect     = [xL xR zB zT];
    out.dx       = dx;
    out.dz       = dz;
end