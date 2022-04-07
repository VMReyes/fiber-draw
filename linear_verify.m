s = tf('s');
H = 1/(s^2+3*s+2);
[MAG,PHASE,W,SDMAG,SDPHASE] = bode(H, {1e-3, 1e2});
t = 0:1e-2:1000; % sampling freq & time simulation length both important
all_w = logspace(-3, 2, 100);
discrete_bode = zeros(size(all_w));
fit_opt = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'Bisquare');
plot_progress = false;

figure(2); latexify_plot;
for i = 1:length(all_w)
    w = all_w(i);
    u = sin(w*t);
    y = lsim(H,u,t);

    eqn_str = sprintf('a*sin(%f*x+phi)+c',w);
    ft = fittype(eqn_str, 'independent', 'x', 'dependent', 'y');
    [fit_obj, goodness_info] = fit(t', y, ft);
    discrete_bode(i) = abs(fit_obj.a);
    
    if plot_progress
        plot(fit_obj, t, y, '-');
    end
    fprintf('%d/%d\n', i, length(all_w))
end

figure(2); latexify_plot; 
loglog(W, MAG(:));
hold on;
loglog(all_w, discrete_bode);
hold off;