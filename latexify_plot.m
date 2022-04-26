set(gcf, 'Color', 'white');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
gr = groot;
for ff = gr.Children
    all_axes = findobj(ff.Children, 'type', 'Axes');
    for ax = all_axes
        try 
            ax.Title.Interpreter = 'latex';
        catch
            continue;
        end
    end
end