% Test plot for truncated distributions
clearvars; clc;
close all; tic;

% Assumptions and notes
% - uses differences in CDFs for Rprobs
% - examine various types of truncation with same mean
% - loop over various mean sampling fractions
% - no sampling from distributions included

% Directory and folder for saving
thisDir = cd; saveFol = 'control batch';

% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Main code - compute variance of controls

% Mean rho (sampling fractions)
rhoMean = [0.9 0.6 0.3 0.1];  nrho = length(rhoMean);
% Mean R and effective R
Rmean = 3; effMean = Rmean*rhoMean;

% Decide number of methods and trajectories
nMeth = 3; 
methNam = {'uniform', 'sub-spread', 'super-spread'};

% Dispersion k of gamma on R
k = logspace(log10(0.5), log10(50), 20); 
lenk = length(k); k10 = log10(k);
% Size of R = length(Iday) + zlen
zlen = 100; % no. zeros to append

% Pdf of R for plotting 
Rpdf = 0.0001:0.0001:50; lenR = length(Rpdf);

% Var to mean ratios and R probabilities
VMs = cell(1, nrho); Rprob0s = VMs; Rprob1s = VMs; Rprob2s = VMs;

for ii = 1:nrho
    % Main code for computing impact of controls on R
    [VM, Rprob0, Rprob1, Rprob2, savStr] = getTruncCtrl(nMeth, k, lenk, Rpdf,...
        effMean(ii), rhoMean(ii), Rmean, saveFol, thisDir, grey1);
    
    % Store loop variables
    VMs{ii} = VM; Rprob0s{ii} = Rprob0;
    Rprob1s{ii} = Rprob1; Rprob2s{ii} = Rprob2;
    disp(['Completed ' num2str(ii) ' of ' num2str(nrho)]);
end

%% Visualisation and processing

% Colours
cols{1} = 'g'; cols{2} = 'b'; cols{3} = 'r';

% Combine all variance to mean ratios
figure;
for ii = 1:nrho
    subplot(ceil(nrho/2), 2, ii);
    plot(k10, VMs{ii}(1,:), 'Color', cols{1}, 'LineWidth', 2);
    hold on;
    plot(k10, VMs{ii}(2,:), 'Color', cols{2}, 'LineWidth', 2);
    plot(k10, VMs{ii}(3,:), 'Color', cols{3}, 'LineWidth', 2);
    hold off;
    xlim([k10(1) k10(end)]);
    h = gca; hL = h.XTickLabel; lenh = length(hL);
    for i = 1:lenh
        h.XTickLabel{i} = round(10^str2double(hL{i}), 2, 'significant');
    end
    xlabel('$k$', 'FontSize', 18);
    ylabel('VM[$R$]', 'FontSize', 18);
    grid off; box off;
    title(['$\rho \mu$ = ' num2str(effMean(ii))], 'FontSize', 18);
    %legend('uniform', 'sub-spread', 'super-spread', 'Location', 'best')
end
cd(saveFol);
saveas(gcf, ['VMcomb_' num2str(Rmean)], 'fig');
cd(thisDir);

% Collect all smallest k distributions
idmin = 1;
P0 = zeros(nrho, lenR); P1 = P0; P2 = P0; 
for ii = 1:nrho
    %P0(ii, :) = cumtrapz(Rpdf, Rprob0s{ii}(2,:));
    P0(ii, :) = cumsum(Rprob0s{ii}(idmin,:));
    %P1(ii, :) = cumtrapz(Rpdf, Rprob1s{ii}(2,:));
    P1(ii, :) = cumsum(Rprob1s{ii}(idmin,:));
    %P2(ii, :) = cumtrapz(Rpdf, Rprob2s{ii}(2,:));
    P2(ii, :) = cumsum(Rprob2s{ii}(idmin,:));
end
R10 = log10(Rpdf); kmin = k(idmin);

% Plot for max and min rho value the 3 controls
figure; jj = 1;
for ii = [1 nrho]
    subplot(1, 2, jj);
    semilogx([10 10], [0 1], '--', 'Color', grey1, 'LineWidth', 2);
    hold on; xlim([0.01 Rpdf(end)])
    semilogx([1/10 1/10], [0 1], '--', 'Color', grey1, 'LineWidth', 2);
    semilogx(Rpdf, P0(ii, :), 'Color', cols{1}, 'LineWidth', 2);
    semilogx(Rpdf, P1(ii, :), 'Color', cols{2}, 'LineWidth', 2);
    semilogx(Rpdf, P2(ii, :), 'Color', cols{3}, 'LineWidth', 2);
    xlabel(['$R\,|\, k = $' num2str(kmin)], 'FontSize', 18);
    ylabel('F($R$)', 'FontSize', 18);
    grid off; box off; hold off;
    title(['$\rho \mu$ = ' num2str(effMean(ii))], 'FontSize', 18);
    jj = jj + 1;
end
cd(saveFol);
saveas(gcf, ['Rcdfcomb_' num2str(Rmean)], 'fig');
cd(thisDir);

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);


