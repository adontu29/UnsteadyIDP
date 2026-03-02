function plotPotentialRegion(filename)
% plotPotentialRegion  Identify and plot "potential flow" region for one file.
%
% Criteria:
%   - isValid == 1
%   - |vorticity| < omega_thr   (threshold from percentile)
%   - |divergence| < div_thr    (threshold from percentile)
%
% Overlay the resulting region on u(x,z).

    % ---------- Load data ----------
    S    = load(filename);
    data = S.datasa;

    % --- Coordinates and velocity (airfoil or lab frame, doesn't matter here) ---
    x = data.v1;    % [mm]
    z = data.v3;    % [mm]
    u = data.v4;    % [m/s]

    % --- Scalars for criteria (adjust v# if needed) ---
    omegaMag = data.v20;   % |Vorticity| [1/s]
    div2D    = data.v21;   % Divergence 2D [1/s]
    isValid  = data.v44;   % 0/1 validity mask

    % ---------- Build "valid" mask ----------
    valid = (isValid == 1) & ...
            ~isnan(omegaMag) & ~isnan(div2D) & ~isnan(u);

    % Absolute values
    absOmega = abs(omegaMag);
    absDiv   = abs(div2D);

    % ---------- Choose thresholds from data ----------
    % Use lower-percentile values (quietest 40%) as "potential" region.
    % Adjust 40 -> 30 or 50 to be stricter/looser.
    omega_thr = prctile(absOmega(valid), 40);
    div_thr   = prctile(absDiv(valid),   80);

    fprintf('File: %s\n', filename);
    fprintf('  omega_thr = %.3g 1/s\n', omega_thr);
    fprintf('  div_thr   = %.3g 1/s\n', div_thr);

    % ---------- Potential-flow mask ----------
    potMask = valid & ...
              (absOmega <= omega_thr) & ...
              (absDiv   <= div_thr);

    % ---------- Plot u(x,z) with potential region outlined ----------
    umin = min(u(valid));
    umax = max(u(valid));
    Nlev = 30;
    levels = linspace(umin, umax, Nlev+1);

    figure;
    contourf(x, z, u, levels);
    clim([umin umax]);
    axis equal tight;
    colormap(parula(Nlev));
    c = colorbar; c.FontSize = 24;
    title(sprintf('u(x,z) and potential region - %s', filename), ...
          'Interpreter','none','FontSize',24);
    xlabel('x [mm]'); ylabel('z [mm]');
    hold on;

    % Draw boundary of potential region (where mask switches from 0 to 1)
    contour(x, z, potMask, [1 1], 'k', 'LineWidth', 2);

    hold off;

    % ---------- Optional: separate mask-only plot ----------
    figure;
    pcolor(x, z, double(potMask));
    shading flat; axis equal tight;
    colormap(gray(2));
    colorbar('Ticks',[0 1],'TickLabels',{'non-potential','potential'});
    title(sprintf('Potential-flow mask - %s', filename), ...
          'Interpreter','none','FontSize',24);
    xlabel('x [mm]'); ylabel('z [mm]');
end