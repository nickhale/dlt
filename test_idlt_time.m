% Test IDLT() with a random vector. Compare with direct approach.

nn = ceil(logspace(3,5,1001));
loopnum = 5;

t_rec(nn(end),1) = 0; t_cheb1(nn(end),1) = 0;
nn = [nn(1), nn];

for n = nn
    fprintf('%i\n', n);   
    c = rand(n, 1);

    % SLOW:
    tic
    for loop = 1:loopnum
        [x_leg, w_leg] = legpts(n);
        scl = ((0:n-1)'+.5);
        wc = (w_leg.'.*c);
%         u = evalP(x_leg)\c;
%         u = scl.*(evalP(x_leg)'*(wc))
        u = scl.*(evalPc_transpose(x_leg, wc));
    end
    t_rec(n) = toc/loopnum;
    
    % FAST:
    tic
    for loop = 1:loopnum
        v = idlt_1(c);
    end
    t_cheb1(n) = toc/loopnum;
end

save test_idlt_time

%%

load test_idlt_time

close all
LW = 'LineWidth';
cols = get(0, 'DefaultAxesColorOrder');

h0 = loglog(NaN, NaN, '-s', NaN, NaN, '-^', LW, 3); hold on
h = loglog(nn, t_rec(nn), 's', nn, t_cheb1(nn), '^', LW, 1);
uistack(h(1), 'top')

nn1000 = nn;
hl = loglog(nn1000, 4.5e-9*nn1000.^2, '--k', ...
    nn1000, 2.7e-7*(nn1000.*log(nn1000).^2./log(log(nn1000))), '--k', LW, 3);

l = legend('recurrence', 'cheb$_1$', 'location', 'NW');
set(l, 'Interpreter', 'LaTeX');

set(gca, 'fontsize', 14)
xlim(nn([1 end]))
ylim([5e-3, 1e1])
grid on

% Choose nicer colors!
% set(h0(3), 'color', cols(5,:));
for k = 1:numel(h)
    set(h(k), 'markerfacecolor', get(h0(k), 'color'));
    set(h(k), 'color', get(h0(k), 'color'));
end
set(h, 'markersize', 1)

nn2 = nn(1:50:end);
hp1 = loglog(nn2, t_cheb1(nn2), '^k', 'linewidth', 1.5);
set(hp1, 'markerfacecolor', get(h(2), 'color'));
hpr = loglog(nn2, t_rec(nn2), 'sk', 'linewidth', 1.5);
set(hpr, 'markerfacecolor', get(h(1), 'color'));

uistack(hp1, 'top')
set(hp1, 'markersize', 6)
set(hpr, 'markersize', 7)

return
%%

pause(1)
print -depsc2 ../paper/figures/idlt_time
