s = tf('s');
H = 1/(s^2+3*s+2); 

plot_progress = false; % true / false
sampling_freq = 100; %Hz, try 10, 100
sim_time = 20; % try 20, 1000
% sampling freq & time simulation length both important

t = 0:1/sampling_freq:sim_time;
[MAG,PHASE,W,SDMAG,SDPHASE] = bode(H, {1e-3, 1e2});
all_w = logspace(-3, 2, 100);
discrete_bode = zeros(size(all_w));
fit_opt = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'Bisquare');


if (plot_progress) figure(1); latexify_plot; end

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
xline(sampling_freq/2)
hold off;
grid on; grid minor;