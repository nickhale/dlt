% Test NDCT() with a random vector. Compare with direct approach.

nn = 1e2:1e4;
loopnum = 5;

t_rec(nn(end),1) = 0; t_chebs(nn(end),1) = 0; t_cheb1(nn(end),1) = 0;
nn = [nn(1), nn];

for n = nn
    fprintf('%i\n', n);
    c = randn(n, 1);
    
%     % SLOW:
%     tic
%     for loop = 1:loopnum
%         x_leg = legpts(n);
%         u = evalT(x_leg)*c;
%     end
%     t_rec(n) = toc/loopnum;
    
    % SLOW:
    tic
    for loop = 1:loopnum
        x_leg = legpts(n);
        u = evalTc(x_leg,c);
    end
    t_rec(n) = toc/loopnum;
        
    % FAST:
    tic
    for loop = 1:loopnum
        v = ndct_1(c);
    end
    t_cheb1(n) = toc/loopnum;
    
    % FAST:
    tic
    for loop = 1:loopnum
        v = ndct_s(c);
    end
    t_chebs(n) = toc/loopnum;

end

save test_ndct_time

%%

load test_ndct_time

close all
LW = 'LineWidth';
cols = get(0, 'DefaultAxesColorOrder');

h0 = loglog(NaN, NaN, '-s', NaN, NaN, '-^', NaN, NaN, '-o', ...
    LW, 3); hold on
h = loglog(nn, t_rec(nn), 's', nn, t_cheb1(nn), '^', nn, t_chebs(nn), 'o', ...
    LW, 1);
uistack(h(1), 'top')

nn2000 = nn(nn>2000);
hl = loglog(nn2000, 4.5e-9*nn2000.^2, '--k', ...
    nn2000, 6.5e-7*(nn2000.*log(nn2000)), '--k', LW, 3);

l = legend('recurrence', 'cheb$_1$', 'cheb$_\ast$', 'location', 'NW');
set(l, 'Interpreter', 'LaTeX');

set(gca,'fontsize',14)
ylim([.3e-3 .4])
% xlim([1e2, 1e4])
grid on, shg

% Choose nicer colors!
set(h0(3), 'color', cols(5,:));
for k = 1:numel(h)
    set(h(k), 'markerfacecolor', get(h0(k), 'color'));
    set(h(k), 'color', get(h0(k), 'color'));
end
set(h, 'markersize', 1)

nn650 = nn(nn>650);
nn1300 = nn(nn>1300);
p1 = polyfit(log(nn650)', log(t_chebs(nn650)), 1);
ps = polyfit(log(nn1300)', log(t_cheb1(nn1300)), 1);
nn2 = logspace(log10(nn(1)), log10(nn(end)), 40);
nn2b = nn2(1:2:end); nn2a = nn2(2:2:end); 
nn2b(nn2b<1300) = []; nn2a(nn2a<650) = [];
hp1 = loglog(nn2a, exp(polyval(p1, log(nn2a))), 'ok', 'linewidth', 1.5);
hps = loglog(nn2b, exp(polyval(ps, log(nn2b))), '^k', 'linewidth', 1.5);
set(hps, 'markerfacecolor', get(h(2), 'color'));
set(hp1, 'markerfacecolor', get(h(3), 'color'));

nn650 = nn(nn<650);
p1 = polyfit(log(nn650)', log(t_chebs(nn650)), 1);
ps = polyfit(log(nn650)', log(t_cheb1(nn650)), 1);
nn2 = ceil(logspace(log10(nn(1)), log10(nn(end)), 40));
nn2b = nn2(1:2:end); nn2a = nn2(2:2:end); 
nn2b(nn2b>1300) = []; nn2a(nn2a>650) = [];
hps = loglog(nn2b, exp(polyval(ps, log(nn2b))), '^k', 'linewidth', 1.5);
hp1 = loglog(nn2a, exp(polyval(p1, log(nn2a))), 'ok', 'linewidth', 1.5);
set(hps, 'markerfacecolor', get(h(2), 'color'));
set(hp1, 'markerfacecolor', get(h(3), 'color'));

hpr = loglog(nn2, t_rec(nn2), 'sk', 'linewidth', 1.5);
set(hpr, 'markerfacecolor', get(h(1), 'color'));

uistack(hps, 'top')
uistack(hp1, 'top')
set(hps, 'markersize', 6)
set(hp1, 'markersize', 6)

return
%%

print -depsc2 ../paper/figures/ndct_time