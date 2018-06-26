function PlotVerticalLines(times, minY, maxY, cmap)

if nargin < 2
    YL = ylim;
    minY = YL(1);
    maxY = YL(2);
end
if nargin < 4
    cmap = 'b-';
end
hold on;
for i= 1:numel(times)
    plot([times(i), times(i)], [minY, maxY], cmap)
end
hold off;