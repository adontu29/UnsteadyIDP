function plotNormal(X,Z,Function,TitleStr,FigureNumber,Start,Step,Stop)

    % --- Handle optional figure number ---
    if nargin < 4 || isempty(FigureNumber)
        figure;      % open a new figure
    else
        figure(FigureNumber);
    end

    NumberOfIntervals = (Stop-Start)/Step;
    contourf(X, Z, Function, Start:Step:Stop)
    clim([Start Stop])
    title(TitleStr, 'FontSize', 32);
    axis equal
    axis tight
    colormap(parula(round(NumberOfIntervals)));
    c = colorbar;
    c.FontSize = 24;

end