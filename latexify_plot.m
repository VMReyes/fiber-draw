set(gcf, 'Color', 'white');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

ff = gcf;
all_axes = findobj(ff.Children, 'type', 'Axes');
for ax = 1:length(all_axes)
    all_axes(ax).FontSize = 10;
    all_axes(ax).Title.Interpreter = 'latex';
    all_axes(ax).Title.FontSize = 12;    
end