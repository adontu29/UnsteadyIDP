function out = solveLaplacePotentialRectExtended(X, Z, u, w, valid, rectExt, rectKeep, Uinf, Winf)
%SOLVELAPLACEPOTENTIALRECTEXTENDED
% Solve Laplace for phi on an upstream-extended rectangular domain.
%
% The measured ROI is rectKeep.
% The solve domain is rectExt, where ONLY the left boundary is moved upstream.
%
% In the artificial upstream strip:
%   - mask = true
%   - velocity is set to uniform freestream (Uinf, Winf)
%
% The solver uses a masked pure-Neumann formulation:
%   dphi/dx = u, dphi/dz = w
% on all outer boundaries / mask boundaries,
% and fixes one reference node phi=0 to remove the nullspace.
%
% Outputs:
%   out.phi_ext      : phi on extended domain
%   out.mask_ext     : mask on extended domain
%   out.Xext, Zext   : extended-domain grids
%   out.phi_keep     : phi only in original measured ROI
%   out.mask_keep    : original measured ROI mask
%   out.keep_idx     : [i1 i2 j1 j2] indices in original measured grid
%   out.keep_rows    : rows in extended domain corresponding to kept ROI
%   out.ref_mask_ext : artificial upstream strip mask for c_dot estimation
%   out.dx, out.dz   : spacings
%
% Assumptions:
%   - X varies mainly along rows
%   - Z varies mainly along columns
%   - rectExt and rectKeep share the same right/top/bottom boundaries
%   - only the left boundary is extended

    arguments
        X (:,:) double
        Z (:,:) double
        u (:,:) double
        w (:,:) double
        valid (:,:) logical
        rectExt (1,4) double
        rectKeep (1,4) double
        Uinf (1,1) double
        Winf (1,1) double
    end

    % ------------------------------------------------------------
    % Unpack rectangles
    % ------------------------------------------------------------
    xL_ext = min(rectExt(1), rectExt(2));
    xR_ext = max(rectExt(1), rectExt(2));
    zB_ext = min(rectExt(3), rectExt(4));
    zT_ext = max(rectExt(3), rectExt(4));

    xL_keep = min(rectKeep(1), rectKeep(2));
    xR_keep = max(rectKeep(1), rectKeep(2));
    zB_keep = min(rectKeep(3), rectKeep(4));
    zT_keep = max(rectKeep(3), rectKeep(4));

    % ------------------------------------------------------------
    % Original grid lines
    % ------------------------------------------------------------
    xline = X(:,1);
    zline = Z(1,:);

    % Find keep-region indices on measured grid
    i1 = find(xline >= xL_keep, 1, 'first');
    i2 = find(xline <= xR_keep, 1, 'last');
    j1 = find(zline >= zB_keep, 1, 'first');
    j2 = find(zline <= zT_keep, 1, 'last');

    if isempty(i1) || isempty(i2) || isempty(j1) || isempty(j2) || i2 < i1 || j2 < j1
        error('rectKeep bounds not found on measured grid.');
    end

    % Check that only left boundary changed
    tol = 5e-12;
    if abs(xR_ext - xR_keep) > tol || abs(zB_ext - zB_keep) > tol || abs(zT_ext - zT_keep) > tol
        error('This solver assumes only the LEFT boundary is extended.');
    end
    if xL_ext > xL_keep + tol
        error('rectExt left boundary must be <= rectKeep left boundary.');
    end

    % ------------------------------------------------------------
    % Extract measured keep region
    % ------------------------------------------------------------
    Xk = X(i1:i2, j1:j2);
    Zk = Z(i1:i2, j1:j2);
    uk = u(i1:i2, j1:j2);
    wk = w(i1:i2, j1:j2);
    mk = valid(i1:i2, j1:j2);

    % Grid spacing
    dx = median(abs(diff(Xk(:,1))));
    dz = median(abs(diff(Zk(1,:))));

    if ~isfinite(dx) || ~isfinite(dz) || dx <= 0 || dz <= 0
        error('Bad dx/dz computed from keep ROI.');
    end

    % ------------------------------------------------------------
    % Build artificial upstream strip
    % ------------------------------------------------------------
    x_keep = Xk(:,1);
    z_keep = Zk(1,:);

    if xL_ext < x_keep(1) - 0.5*dx
        x_art = (x_keep(1)-dx : -dx : xL_ext);
        x_art = fliplr(x_art).';
    else
        x_art = zeros(0,1);
    end

    nArt = numel(x_art);

    x_ext = [x_art; x_keep];
    z_ext = z_keep;

    [Xext, Zext] = ndgrid(x_ext, z_ext);

    nr = size(Xext,1);
    nc = size(Xext,2);

    % ------------------------------------------------------------
    % Assemble extended fields
    % ------------------------------------------------------------
    uext = Uinf * ones(nr, nc);
    wext = Winf * ones(nr, nc);
    mext = true(nr, nc);

    % Insert measured keep-region data in the right part
    uext(nArt+1:end, :) = uk;
    wext(nArt+1:end, :) = wk;
    mext(nArt+1:end, :) = mk;

    % Reference region for c_dot estimation = artificial strip
    ref_mask_ext = false(nr, nc);
    if nArt > 0
        ref_mask_ext(1:nArt, :) = true;
    else
        % fallback if no extension was actually added:
        % use leftmost 5% of kept domain
        nc_ref = max(2, round(0.05 * nc));
        ref_mask_ext(nArt+1:end, 1:nc_ref) = true;
    end
    ref_mask_ext = ref_mask_ext & mext;

    % ------------------------------------------------------------
    % Unknown numbering
    % ------------------------------------------------------------
    id = zeros(nr, nc);
    lin = find(mext);
    N = numel(lin);

    if N < 10
        error('Too few valid points in extended domain (%d).', N);
    end

    id(lin) = 1:N;

    A = spalloc(N, N, 5*N);
    b = zeros(N,1);

    idx = @(ii,jj) id(ii,jj);

    % ------------------------------------------------------------
    % Build masked Laplace system with Neumann BC from local velocity
    %
    % (phiE - 2phiP + phiW)/dx^2 + (phiN - 2phiP + phiS)/dz^2 = 0
    %
    % Missing neighbor => local Neumann BC using ghost-point elimination.
    % ------------------------------------------------------------
    for ii = 1:nr
        for jj = 1:nc
            p = idx(ii,jj);
            if p == 0
                continue;
            end

            cP = -2*(1/dx^2 + 1/dz^2);

            % West neighbor (ii-1, jj)
            if ii-1 >= 1 && idx(ii-1,jj) ~= 0
                q = idx(ii-1,jj);
                A(p,q) = A(p,q) + 1/dx^2;
            else
                b(p) = b(p) - uext(ii,jj)/dx;
                cP = cP + 1/dx^2;
            end

            % East neighbor (ii+1, jj)
            if ii+1 <= nr && idx(ii+1,jj) ~= 0
                q = idx(ii+1,jj);
                A(p,q) = A(p,q) + 1/dx^2;
            else
                b(p) = b(p) + uext(ii,jj)/dx;
                cP = cP + 1/dx^2;
            end

            % South neighbor (ii, jj-1)
            if jj-1 >= 1 && idx(ii,jj-1) ~= 0
                q = idx(ii,jj-1);
                A(p,q) = A(p,q) + 1/dz^2;
            else
                b(p) = b(p) - wext(ii,jj)/dz;
                cP = cP + 1/dz^2;
            end

            % North neighbor (ii, jj+1)
            if jj+1 <= nc && idx(ii,jj+1) ~= 0
                q = idx(ii,jj+1);
                A(p,q) = A(p,q) + 1/dz^2;
            else
                b(p) = b(p) + wext(ii,jj)/dz;
                cP = cP + 1/dz^2;
            end

            A(p,p) = A(p,p) + cP;
        end
    end

    % ------------------------------------------------------------
    % Fix nullspace of pure Neumann problem
    % ------------------------------------------------------------
    ref = 1;
    A(ref,:) = 0;
    A(ref,ref) = 1;
    b(ref) = 0;

    % ------------------------------------------------------------
    % Solve
    % ------------------------------------------------------------
    phi_vec = A \ b;

    phi_ext = nan(nr, nc);
    phi_ext(lin) = phi_vec;

    % ------------------------------------------------------------
    % Keep only original measured ROI part
    % ------------------------------------------------------------
    keep_rows = (nArt+1):nr;

    phi_keep = phi_ext(keep_rows, :);
    mask_keep = mk;

    % only retain measured valid region
    phi_keep(~mask_keep) = NaN;

    % ------------------------------------------------------------
    % Pack output
    % ------------------------------------------------------------
    out = struct();
    out.phi_ext      = phi_ext;
    out.mask_ext     = mext;
    out.Xext         = Xext;
    out.Zext         = Zext;

    out.phi_keep     = phi_keep;
    out.mask_keep    = mask_keep;

    out.keep_idx     = [i1 i2 j1 j2];
    out.keep_rows    = keep_rows;

    out.ref_mask_ext = ref_mask_ext;

    out.dx = dx;
    out.dz = dz;
    out.nArt = nArt;
end