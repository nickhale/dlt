% clc, close all
% 
% nn = 100:1000;
% 
% cols = get(0, 'DefaultAxesColorOrder');
% ell = 0;
% 
% for dec = [0 .5 1 1.5]
%     
%     for n = nn
%         
%         disp(n)
%         rng(0, 'twister');
%         r = randn(n, 1);
%         if ( dec == 0 )
%             c = r;
%         elseif ( dec == .5 )
%             c = r./sqrt(1:n).';
%         elseif ( dec == 1 )
%             c = r./(1:n).';
%         elseif ( dec == 1.5 )
%             c = r./(1:n).'.^1.5;
%         else
%             error
%         end
%         
%         %% SLOW:
%         [x, ~, ~, t] = legpts(n);
%         u = evalTc(x,c);
%         
%         %% QUAD
%         v = ndct_quad(c);
%         
%         %% FAST:
%         w = ndct_s(c);
%         w1 = ndct_1(c);
%         
%         %% TESTING:
%         err_rec(n) = norm((u - v), inf);
%         err_new(n) = norm((w - v), inf);
%         err_new1(n) = norm((w1 - v), inf);
%         
%     end
%     
%     clear ell
%     if ( dec == 0 )
%         save ../paper/figures/ndct_nodecay
%     elseif ( dec == 0.5 )
%         save ../paper/figures/ndct_05decay
%     elseif ( dec == 1 )
%         save ../paper/figures/ndct_10decay
%     elseif ( dec == 1.5 )
%         save ../paper/figures/ndct_15decay
%     end
%     
% end
% 
% return

%%

close all
ell = 0;
for dec = [0 .5 1 1.5]
    
    if ( dec == 0 )
        load ../paper/figures/ndct_nodecay
    elseif ( dec == 0.5 )
        load ../paper/figures/ndct_05decay
    elseif ( dec == 1 )
        load ../paper/figures/ndct_10decay
    elseif ( dec == 1.5 )
        load ../paper/figures/ndct_15decay
    end

    ell = ell + 1;    
    LW = 'LineWidth'; MS = 'MarkerSize';
    for dec2 = [0 .5 1 1.5]
        for f = 1:3
            figure(f)
            loglog([10 100], [1e-14 1e-7], '-', LW, 3); hold on
        end
    end
    
    figure(1)
    h1 = loglog(nn, err_rec(nn), '.', LW, 3);
    set(h1, 'color', cols(ell,:));
    set(h1, 'MarkerFacecolor', get(h1, 'color'));
    figure(2)
    h2 = loglog(nn, err_new(nn), '.', LW, 3);
    set(h2, 'color', cols(ell,:));
    set(h2, 'MarkerFacecolor', get(h2, 'color'));
    figure(3)
    h2b = loglog(nn, err_new1(nn), '.', LW, 3); 
    set(h2b, 'color', cols(ell,:));
    set(h2b, 'MarkerFacecolor', get(h2b, 'color'));
    set(h1, MS, 12), set(h2, MS, 12),set(h2b, MS, 12)

    for f = 1:3
        figure(f), xlim([100, 1000]), ylim([1e-16 1e-9]), grid on
    end
    if ( dec == 0 )
        figure(1), h3 = loglog(nn, 2e-10*(nn/nn(end)).^2.5, '--k', LW, 3);
        figure(2), h4 = loglog(nn, 5e-16*nn.^2./log(nn).^2, '--k', LW, 3);
        figure(3), h5 = loglog(nn, 5e-16*nn.^2./log(nn).^2, '--k', LW, 3);
    elseif ( dec == 0.5 )
        figure(1), h3 = loglog(nn, 8e-18*(nn).^2, '--k', LW, 3);
        figure(2), h4 = loglog(nn, 6.3e-16*nn.^1.5./log(nn).^2, '--k', LW, 3);
        figure(3), h5 = loglog(nn, 6.3e-16*nn.^1.5./log(nn).^2, '--k', LW, 3);
    elseif ( dec == 1 )
        figure(1), h3 = loglog(nn, 8.4e-18*nn.^1.5, '--k', LW, 3);
        figure(2), h4 = loglog(nn, 8.3e-17*sqrt(nn).*log(nn), '--k', LW, 3);
        figure(3), h5 = loglog(nn, 8.3e-17*sqrt(nn).*log(nn), '--k', LW, 3);
    elseif ( dec == 1.5 )
        figure(1), h3 = loglog(nn, 1.4e-17*nn, '--k', LW, 3);
        figure(2), h4 = loglog(nn, 5.5e-17*log(nn).^2, '--k', LW, 3);
        figure(3), h5 = loglog(nn, 5.5e-17*log(nn).^2, '--k', LW, 3);
    end

    HV = 'HandleVisibility';
    set(h1, HV, 'off'), set(h2, HV, 'off'), set(h2b, HV, 'off'), 
    set(h3, HV, 'off'), set(h4, HV, 'off'), set(h5, HV, 'off')
    
end

for f = 1:3
    figure(f)
    l = legend('$\mathcal{O}(n^0)$', '$\mathcal{O}(n^{-0.5})$', ...
        '$\mathcal{O}(n^{-1})$', '$\mathcal{O}(n^{-1.5})$', 'location', 'NW');
    set(l, 'Interpreter', 'LaTeX');
    set(gca, 'fontsize', 14)
end

return
%%

pause(1)
figure(1)
print -depsc2 ../paper/figures/ndct_all_rec
figure(2)
print -depsc2 ../paper/figures/ndct_all_new
figure(2)
print -depsc2 ../paper/figures/ndct_all_dir
