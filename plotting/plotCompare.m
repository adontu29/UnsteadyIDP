function plotCompare(X,Z,Function1,Function2,Title1,Title2,FigureNumber,Start,Step,Stop)

if isempty(FigureNumber)
    f = figure;      % open a new figure
else
    f = figure(FigureNumber);clf;
end


f.WindowState = 'maximized';
f.InnerPosition = [1000,1000,1000,1000];
NumberOfIntervals = (Stop-Start)/Step;
subplot(2,1,1)
contourf(X,Z,Function1,[Start:Step:Stop])
clim([Start Stop])
title(Title1);
axis equal
axis tight
colormap(jet(round(NumberOfIntervals)));
colorbar

subplot(2,1,2)
contourf(X,Z,Function2,[Start:Step:Stop])
clim([Start Stop])
title(Title2);
axis equal
axis tight
colormap(jet(round(NumberOfIntervals)));
colorbar
end