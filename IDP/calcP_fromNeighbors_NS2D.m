function LocalPressure = calcP_fromNeighbors_NS2D(X, Z, U, W, p, Mask, PX, PZ, i, k, weightingMode)

LocalPressure = NaN;

SumWeightedP = 0;
SumWeights   = 0;

neighbors = [
    -1,  0;
     1,  0;
     0,  1;
     0, -1;
    -1, -1;
    -1,  1;
     1, -1;
     1,  1
];

useIso = strcmp(weightingMode,'isopotential');
if useIso
    t0 = [-W(i,k), U(i,k)];   % tangent direction (same idea as your old)
    nt0 = norm(t0);
    if nt0 > 0
        t0 = t0/nt0;
    else
        useIso = false;
    end
end

[sizeX,sizeZ] = size(X);

for n = 1:size(neighbors,1)
    ni = i + neighbors(n,1);
    nk = k + neighbors(n,2);

    if ni < 1 || ni > sizeX || nk < 1 || nk > sizeZ
        continue;
    end
    if Mask(ni,nk) ~= 1
        continue;
    end
    if ~isfinite(p(ni,nk))
        continue;
    end

    dx = (X(ni,nk) - X(i,k));
    dz = (Z(ni,nk) - Z(i,k));
    dist = hypot(dx,dz);
    if dist == 0
        continue;
    end

    unitDir = [dx dz] / dist;

    gradP_cent = [PX(i,k),   PZ(i,k)];
    gradP_nghb = [PX(ni,nk), PZ(ni,nk)];
    gradP = 0.5*(gradP_cent + gradP_nghb);

    if any(~isfinite(gradP))
        continue;
    end

    dP = dot(gradP, unitDir);

    AdjustedPressure = p(ni,nk) - dP*dist;

    if useIso
        dx_abs = abs(X(ni,nk) - X(i,k));
        dz_abs = abs(Z(ni,nk) - Z(i,k));
        dist_abs = hypot(dx_abs, dz_abs);
        if dist_abs == 0, continue; end
        stepVec = [dx_abs dz_abs]/dist_abs;
        wgt = abs(dot(stepVec, t0))^2;
    else
        wgt = 1.0;
    end

    SumWeightedP = SumWeightedP + wgt*AdjustedPressure;
    SumWeights   = SumWeights   + wgt;
end

if SumWeights > 0
    LocalPressure = SumWeightedP / SumWeights;
end

end