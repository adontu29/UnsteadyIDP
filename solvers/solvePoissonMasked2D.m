function phi = solvePoissonMasked2D(X, Z, f, potMask, phiFixed)

[Nz, Nx] = size(X);
dx = abs(X(2,1) - X(1,1))/1000;    
dz = abs(Z(1,2) - Z(1,1))/1000;   

N = Nx * Nz;
idx = @(i,j) sub2ind([Nz Nx], i, j);

A = spalloc(N, N, 5*N);
b = zeros(N,1);

for j = 1:Nx
    for i = 1:Nz

        k = idx(i,j);

        % If node is outside potential region 
        if ~potMask(i,j)
            A(k,k) = 1;
            b(k)   = phiFixed;
            continue;
        end

        % Outer physical boundary: also treat as Dirichlet
        if i==1 || i==Nz || j==1 || j==Nx
            A(k,k) = 1;
            b(k)   = phiFixed;
            continue;
        end

        % Interior node in Ω_p: 5-pt stencil, but neighbors
        % outside Ω_p treated as Dirichlet = phiFixed.
        % Start with center coefficient:
        coeff_center = -2*(1/dx^2 + 1/dz^2);
        rhs_k = f(i,j);

        % Neighbor indices
        neigh = [i, j+1;
                 i, j-1;
                 i+1, j;
                 i-1, j];

        for n = 1:4
            ii = neigh(n,1);
            jj = neigh(n,2);
            kk = idx(ii,jj);
            if potMask(ii,jj)
                % neighbor in Ω_p -> regular stencil
                if n <= 2   % left/right
                    coeff = 1/dx^2;
                else        % up/down
                    coeff = 1/dz^2;
                end
                A(k,kk) = coeff;
            else
                % neighbor outside Ω_p -> Dirichlet to phiFixed
                if n <= 2
                    coeff = 1/dx^2;
                else
                    coeff = 1/dz^2;
                end
                % move known neighbor to RHS: coeff * phiFixed
                rhs_k = rhs_k - coeff * phiFixed;
            end
        end

        A(k,k) = coeff_center;
        b(k)   = rhs_k;

    end
end

phi_vec = A \ b;
phi = reshape(phi_vec, Nz, Nx);

end