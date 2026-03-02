function plotVelSlices(filename)
% plotVelSlices  Plot u-velocity field and airfoil outline for one dataset
%                using a style similar to plotNormal.

    % Load file and get inner struct
    S     = load(filename);
    data  = S.datasa;      % everything is inside 'datasa'

    ['AOA is ' num2str(data.AOA) ' and theta is ' num2str(data.theta)]

    % --- Flow field variables (adapt if your mapping is different)
    x = data.v45;    % [mm]
    y = data.v2;    % [mm]
    z = data.v46;    % [mm]
    u = data.v4;    % [m/s]
    v = data.v5;    % [m/s]
    w = data.v6;    % [m/s]

    VelMag = sqrt(u.^2+w.^2);

    % --- Airfoil coordinates (2 x N: first row x, second row z)
    x_airfoil = data.airfoilrot(1, :);
    z_airfoil = data.airfoilrot(2, :);

    % Number of intervals (you can tweak this)
    NumberOfIntervals = 30;
    
    Start = 0;
    Stop  = 20;
    Step  = (Stop - Start) / NumberOfIntervals;

    % --- Plot u(x,z) with contourf in the same style ---
    figure;
    contourf(x, z, VelMag, Start:Step:Stop);
    clim([Start Stop]);
    axis equal;
    axis tight;
    colormap(parula(round(NumberOfIntervals)));
    c = colorbar;
    c.FontSize = 24;
    title(['AOA' num2str(data.AOA)], 'Interpreter', 'none', 'FontSize', 32);
    xlabel('x [mm]');
    ylabel('z [mm]');
    hold on;

    % --- Overlay airfoil outline ---
    if ~isempty(x_airfoil) && ~isempty(z_airfoil)
        plot(x_airfoil, z_airfoil, 'k-', 'LineWidth', 2);
    end

    hold off;

end