ccc

nn = 10:1:1000;

cols = get(0, 'DefaultAxesColorOrder');

for n = nn
    
    disp(n)
    rng(0, 'twister');
    c = randn(n, 1);
    
    %% SUPER SLOW:
%     [x, ~, ~, t] = legpts(n);
%     u = evalP(x)\c;
    
    %% SLOW
    [x, w, ~, t] = legpts(n);
    ub = evalPc_transpose(x, (w'.*c)) .* ((0:n-1)'+.5);
    
    %% QUAD
    v = idlt_quad(c);
    
    %% FAST:
    w = idlt_1(c);
    
    %% TESTING:
%     err_rec(n) = norm((u - v), inf);
    err_rec(n) = norm((ub - v), inf);
    err_new(n) = norm((w - v), inf);
    
end

save ../paper/figures/idlt_nodecayb
% return
%%
close all
load ../paper/figures/idlt_nodecay


h0 = loglog(NaN, NaN, '.-', 'LineWidth', 3); hold on
set(h0, 'color', cols(1,:));
h0 = loglog(NaN, NaN, '.-', 'LineWidth', 3); hold on
set(h0, 'color', cols(2,:));

h2 = loglog(nn, err_new(nn), '.', 'LineWidth', 3);
set(h2, 'color', cols(2,:));
set(h2, 'MarkerFacecolor', get(h2, 'color'));

h1 = loglog(nn, err_rec(nn), '.', 'LineWidth', 3); hold on
set(h1, 'color', cols(1,:));
set(h1, 'MarkerFacecolor', get(h1, 'color'));

% h3 = loglog(nn, 6e-16*(nn), '--k', 'LineWidth', 3); hold on
h3 = loglog(nn, 1e-16*nn.*log(nn), '--k', 'LineWidth', 3); hold on
set(h3, 'HandleVisibility', 'off')

set(h1, 'Markersize', 12)
set(h2, 'Markersize', 12)

xlim([100, 1000])
ylim([1e-14 3e-12])

grid on
l = legend('recurrence', 'cheb$_1$', 'location', 'NW');
set(l, 'Interpreter', 'LaTeX');
set(gca, 'fontsize', 14), shg
return

%%

pause(1)

print -depsc2 ../paper/figures/idlt_results

alignfigs
