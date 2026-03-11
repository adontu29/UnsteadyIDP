function plotVelSlices(filename, mode)
% plotVelSlices  Plot u-velocity field and airfoil outline for one dataset
%                using a style similar to plotNormal.

% Load file and get inner struct
S     = load(filename);
data  = S.datasa;      % everything is inside 'datasa'

['AOA is ' num2str(data.AOA) ' and theta is ' num2str(data.theta)]

% --- Flow field variables (adapt if your mapping is different)
x = flip(data.v1,2);    % [mm]
z = flip(data.v3,2);    % [mm]
u = flip(data.v4,2);
w = -flip(data.v6,2);

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
if mode=="contour"
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

else
    figure
    quiver(x,z,u,w, 2)
    hold on
end

% 
% % --- Overlay airfoil outline ---
% if ~isempty(x_airfoil) && ~isempty(z_airfoil)
%     plot(x_airfoil, z_airfoil, 'k-', 'LineWidth', 2);
% end


hold off;

end