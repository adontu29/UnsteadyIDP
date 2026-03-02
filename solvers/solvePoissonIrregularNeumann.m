function phi = solvePoissonIrregularNeumann(X, Z, f, validMask, u, w, phiFixed)
% solvePoissonIrregularNeumann
%
% Solve 2D Poisson / Laplace equation on an irregular domain defined by
% validMask, with Neumann boundary conditions derived from measured
% velocities u,w on the boundary, and a Dirichlet "gauge" in each
% connected component to remove the constant nullspace.
%
%   ∇² phi = f      in {validMask == 1}
%   ∂phi/∂n = u·n   on boundary (approximated by one-sided FD)
%
% Inputs:
%   X,Z        : (Nz x Nx) grid coordinates (regular grid)
%   f          : RHS (same size) – use 0 for Laplace, or divergence for Poisson
%   validMask  : logical (Nz x Nx), 1 = inside measurement domain
%   u,w        : velocity components on same grid (for Neumann BCs)
%   phiFixed   : scalar, reference potential value (gauge)
%
% Output:
%   phi        : potential (Nz x Nx); only meaningful where validMask == 1

    [Nz, Nx] = size(X);

    % Grid spacing (must match main script!)
    dx = abs(X(2,1) - X(1,1))/1000;
    dz = abs(Z(1,2) - Z(1,1))/1000;

    N   = Nx * Nz;
    idx = @(i,j) sub2ind([Nz Nx], i, j);

    A = spalloc(N, N, 5*N);
    b = zeros(N,1);

    for j = 1:Nx
        for i = 1:Nz

            k = idx(i,j);

            % ---- Outside measurement domain: dummy Dirichlet node ----
            if ~validMask(i,j)
                A(k,k) = 1;
                b(k)   = phiFixed;
                continue;
            end

            % Check 4-neighbour validity (treat outside array as invalid)
            left_in  = (j > 1   && validMask(i,j-1));
            right_in = (j < Nx  && validMask(i,j+1));
            down_in  = (i > 1   && validMask(i-1,j));  % lower index
            up_in    = (i < Nz  && validMask(i+1,j));  % higher index

            isInterior = left_in && right_in && up_in && down_in;

            % -----------------------------------------------------------------
            % Boundary node: at least one neighbour is outside "validMask"
            % -----------------------------------------------------------------
            if ~isInterior

                % Prefer x-boundary if present, otherwise z-boundary.
                if ~left_in && right_in
                    % Open to the left, normal approx n = (-1,0)
                    % (phi(i,j+1) - phi(i,j))/dx = -u(i,j)
                    kk = idx(i, j+1);
                    A(k,k)  = -1/dx;
                    A(k,kk) =  1/dx;
                    b(k)    = -u(i,j);

                elseif ~right_in && left_in
                    % Open to the right, n ≈ (1,0)
                    % (phi(i,j) - phi(i,j-1))/dx = u(i,j)
                    kk = idx(i, j-1);
                    A(k,k)  =  1/dx;
                    A(k,kk) = -1/dx;
                    b(k)    =  u(i,j);

                elseif ~down_in && up_in
                    % Open "downwards", n ≈ (0,-1)
                    % (phi(i+1,j) - phi(i,j))/dz = -w(i,j)
                    kk = idx(i+1, j);
                    A(k,k)  = -1/dz;
                    A(k,kk) =  1/dz;
                    b(k)    = -w(i,j);

                elseif ~up_in && down_in
                    % Open "upwards", n ≈ (0,1)
                    % (phi(i,j) - phi(i-1,j))/dz = w(i,j)
                    kk = idx(i-1, j);
                    A(k,k)  =  1/dz;
                    A(k,kk) = -1/dz;
                    b(k)    =  w(i,j);

                else
                    % Pathological case: isolated pixel or weird corner.
                    % Just pin to gauge.
                    A(k,k) = 1;
                    b(k)   = phiFixed;
                end

                continue;  % done with this boundary node
            end

            % -----------------------------------------------------------------
            % Interior node: full 5-point Poisson stencil
            % -----------------------------------------------------------------
            coeff_center = -2*(1/dx^2 + 1/dz^2);
            rhs_k        = f(i,j);

            % right neighbour
            kk = idx(i, j+1);
            A(k,kk) = 1/dx^2;

            % left neighbour
            kk = idx(i, j-1);
            A(k,kk) = 1/dx^2;

            % up neighbour
            kk = idx(i+1, j);
            A(k,kk) = 1/dz^2;

            % down neighbour
            kk = idx(i-1, j);
            A(k,kk) = 1/dz^2;

            A(k,k) = coeff_center;
            b(k)   = rhs_k;

        end
    end

       % -----------------------------------------------------------------
    % Gauge fix: one node per connected component in validMask
    % -----------------------------------------------------------------
    labels = bwlabel(validMask, 4);       % 4-connectivity
    nComp  = max(labels(:));

    for c = 1:nComp
        compIdx = find(labels == c, 1, 'first');   % linear index in matrix sense
        if ~isempty(compIdx)
            k0 = compIdx;    % same linear indexing as A's rows
            A(k0,:) = 0;
            A(k0,k0) = 1;
            b(k0)   = phiFixed;
        end
    end

    % -----------------------------------------------------------------
    % Strong diagonal regularisation to kill remaining near-nullspace
    % -----------------------------------------------------------------
    diagA   = abs(diag(A));
    maxDiag = max(diagA);
    if maxDiag == 0
        maxDiag = 1;   % safety
    end
    lambda = 1e-3 * maxDiag;   % you can tune 1e-3 → 1e-2 if needed

    A = A + lambda * speye(N);

    % -----------------------------------------------------------------
    % Solve linear system
    % -----------------------------------------------------------------
    phi_vec = A \ b;
    phi     = reshape(phi_vec, Nz, Nx);

    % Optional: remove mean inside valid domain (pure gauge)
    mphi = mean(phi(validMask), 'omitnan');
    phi  = phi - mphi;
end