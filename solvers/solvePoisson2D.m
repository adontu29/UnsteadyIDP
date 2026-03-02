function phi = solvePoisson2D(X, Z, f, phiBC)

[Nz, Nx] = size(X);

% Grid spacing (assume uniform)
dx = mean(diff(X(1,:)));
dz = mean(diff(Z(:,1)));

% Number of unknowns
N = Nx * Nz;

% map (i,j) to linear index
idx = @(i,j) sub2ind([Nz Nx], i, j);

% Sparse matrix and RHS vector
A = spalloc(N, N, 5*N);
b = zeros(N,1);

for j = 1:Nx
    for i = 1:Nz

        k = idx(i,j);

        if i==1 || i==Nz || j==1 || j==Nx
            % Boundary point Dirichlet
            A(k,k) = 1;
            b(k)   = phiBC(i,j);

        else
            % Interior point: 5-point Laplacian
            A(k,idx(i,j))     = -2*(1/dx^2 + 1/dz^2);
            A(k,idx(i,j+1))   =  1/dx^2;
            A(k,idx(i,j-1))   =  1/dx^2;
            A(k,idx(i+1,j))   =  1/dz^2;
            A(k,idx(i-1,j))   =  1/dz^2;

            b(k) = f(i,j);
        end
    end
end

% Solve linear system
phi_vec = A \ b;

% Reshape back to 2D
phi = reshape(phi_vec, Nz, Nx);

end