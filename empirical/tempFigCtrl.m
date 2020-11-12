% Reconstruct Rpdf figure temporarily 

% Load largest R then do this
load('/Users/kp10/Desktop/Imperial/2019/Term 3/Code/Super spreading end/sub_spread_code/variance/control batch/ctrlVar_2x7_3_20_3.mat')
ii = 1;
% Other variables
lenR = length(Rpdf); nrho = 2; grey2 = 0.5*ones(1, 3);
P0s = zeros(nrho, lenR); P1s = P0s; P2s = P0s;
VMs = cell(1, nrho); Rprob0s = VMs; Rprob1s = VMs; Rprob2s = VMs;
effMeans = zeros(1, nrho); effMeans(ii) = effMean;

Rprob0s{ii} = Rprob0;
Rprob1s{ii} = Rprob1; 
Rprob2s{ii} = Rprob2;

Rprob1s{ii} = Rprob1; Rprob2s{ii} = Rprob2;
P0s(ii, :) = cumsum(Rprob0s{ii}(2,:));
P1s(ii, :) = cumsum(Rprob1s{ii}(2,:));
P2s(ii, :) = cumsum(Rprob2s{ii}(2,:));

% Then load smallest R
load('/Users/kp10/Desktop/Imperial/2019/Term 3/Code/Super spreading end/sub_spread_code/variance/control batch/ctrlVar_0x3_3_20_3.mat');
ii = 2;
Rprob0s{ii} = Rprob0;
Rprob1s{ii} = Rprob1; 
Rprob2s{ii} = Rprob2;

Rprob1s{ii} = Rprob1; Rprob2s{ii} = Rprob2;
P0s(ii, :) = cumsum(Rprob0s{ii}(2,:));
P1s(ii, :) = cumsum(Rprob1s{ii}(2,:));
P2s(ii, :) = cumsum(Rprob2s{ii}(2,:));
effMeans(ii) = effMean;


% Plot results
figure; jj = 1;
for ii = [1 2]
    subplot(1, 2, jj);
    semilogx(Rpdf, P0s(ii, :), 'LineWidth', 2);
    hold on; xlim([0.05 Rpdf(end)])
    plot(Rpdf, P1s(ii, :), 'LineWidth', 2);
    plot(Rpdf, P2s(ii, :), 'LineWidth', 2);
    plot([10 10], [0 1], '--', 'Color', grey2, 'LineWidth', 2);
    plot([1/10 1/10], [0 1], '--', 'Color', grey2, 'LineWidth', 2);
    xlabel('$R$', 'FontSize', 18);
    ylabel('F($R$)', 'FontSize', 18);
    grid off; box off; hold off;
    title(['$\rho \mu$ = ' num2str(effMeans(ii))], 'FontSize', 18);
    jj = jj + 1; ylim([0 1.01]);
end
