% Batch across elimination probabilities
clearvars; clc;
close all; tic;

% Assumptions and notes
% - reduced R is what expect to affect future incidence
% - treat control and undersampling in same sense
% - fix Rmean to investigate effect of control via rhoMean
% - either pop wide or indiv based with upper or lower truncation
% - R sampled from a distribution over time and trajectories
% - declarations are relative to tst; when last case was seen
% - incidence curve from R (Ebola) forms the data

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));
% Mex code for faster Lam computations
addpath(genpath('/Users/kp10/Desktop/Imperial/2019/Term 3/Code/Super-spreading/super_end/codegen/mex'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', 18);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Save figures
saveFig = 0;
% Folders for loading/saving
loadFol = 'mers'; saveFol = 'mers_batch';

%% Main inputs to elimVaryRVaryk and results

% Confidence level for declaration
mu = 0.95; disp(['Confidence = ' num2str(mu)]);
% Mean R values to examine
Rmean = 1; disp(['Mean R = ' num2str(Rmean)]);

% Three methods of control/undersampling
ctrlNames = {'constant', 'lower', 'upper'};
ctrlMeth = 1; disp(['Control/sampling method = ' ctrlNames(ctrlMeth)]);

% Dispersion k of gamma on R
k = logspace(-1, 1, 15); lenk = length(k);
% Sample fraction
rhoMean = linspace(0.5, 1, 4); lenrho = length(rhoMean);  

% Output cells
z0Mean = cell(1, 1); z0High = z0Mean; z0Low = z0Mean;
z0Worst = z0Mean; tdecAll = z0Mean;

% Main results
parfor j = 1:lenrho
    % Prob of elimination and declaration stats
    [z0Mean{j}, z0High{j}, z0Low{j}, z0Worst{j}, tdecAll{j}] = elimProbCtrl(k, Rmean, mu, rhoMean(j),...
        ctrlMeth, [saveFol num2str(ctrlMeth)], loadFol);
    disp('////////////////////////////////////////////');
    disp(['Completed = ' num2str(j) ' of ' num2str(lenrho)]);
    disp('////////////////////////////////////////////');
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Current time of simulation
simDate = string(datetime('now'));

% Save batch data
save(join(['batchCtrl_' num2str(ctrlMeth) '_' num2str(lenk) '_' num2str(lenrho) '_' simDate '.mat']));


%% Processing results per Rmean

% For plotting
kplt = round(k, 2, 'significant');
Rmeanplt = round(Rmean, 2, 'significant');
rhoplt = round(rhoMean, 2, 'significant');
% Colours for lenrho
cols = cell(1, lenrho); 
cols(cellfun(@isempty, cols)) = {grey1};
cols{1} = 'r'; cols{end} = 'b';
% Colours for lenk
cols2 = cell(1, lenk);
cols2(cellfun(@isempty, cols2)) = {grey1};
cols2{1} = 'r'; cols2{end} = 'b';

% Collect mean and worst tdec at every k
tmeank = cell(1, lenk); tworstk = tmeank;
for i = 1:lenk
    ttemp = zeros(1, lenrho); ttemp2 = ttemp;
    for j = 1:lenrho
        ttemp(j) = tdecAll{j}(i, 1);
        ttemp2(j) = tdecAll{j}(i, end);
    end
    tmeank{i} = ttemp; tworstk{i} = ttemp2;
end
tmeank = cell2mat(tmeank'); tworstk = cell2mat(tworstk');

% Examine trends against k for each rho
figure;
semilogx(k, tmeank(:, 2:end-1), 'Color', grey1, 'LineWidth', 2);
hold on;
semilogx(k, tmeank(:, 1), 'Color', 'r', 'LineWidth', 2);
semilogx(k, tmeank(:, end), 'Color', 'b', 'LineWidth', 2);
semilogx(k, tworstk(:, 2:end-1), 'Color', grey1, 'LineWidth', 2);
semilogx(k, tworstk(:, 1), 'Color', 'r', 'LineWidth', 2);
semilogx(k, tworstk(:, end), 'Color', 'b', 'LineWidth', 2);
xlabel('$k$'); ylabel('$t_{95}$');
grid off; hold off; box off;
if saveFig
    saveas(gcf, ['tdecOverk_' num2str(lenk) '_' num2str(lenrho)], 'fig');
end

% All mean z0 curves with rho and k for given R
figure;
for i = 1:lenk
    subplot(ceil(lenk/2), 2, i);
    hold on;
    for j = 1:lenrho
        plot(z0Mean{j}(i, :), 'Color', cols{j}, 'LineWidth', 2);
    end
    hold off; grid off; box off;
    ylabel(['$\bar{z}_s \, | \, k = $' num2str(kplt(i))]);
    if any(i == [lenk-1, lenk])
        xlabel('$s$ (days)');
    end
end
if saveFig
    saveas(gcf, ['allz0MeanRho_' num2str(lenk) '_' num2str(lenrho)], 'fig');
end

% All mean z0 curves with rho and k for given R
figure;
for i = 1:lenk
    subplot(ceil(lenk/2), 2, i);
    hold on;
    for j = 1:lenrho
        plot(z0Worst{j}(i, :), 'Color', cols2{j}, 'LineWidth', 2);
    end
    hold off; grid off; box off;
    ylabel(['$z_{s, min} \, | \, k = $' num2str(kplt(i))]);
    if any(i == [lenk-1, lenk])
        xlabel('$s$ (days)');
    end
end
if saveFig
    saveas(gcf, ['allz0WorstRho_' num2str(lenk) '_' num2str(lenrho)], 'fig');
end

% All mean z0 curves with k for various rho
figure;
for i = 1:lenrho
    subplot(ceil(lenrho/2), 2, i);
    hold on;
    for j = 1:lenk
        plot(z0Mean{i}(j, :), 'Color', cols2{j}, 'LineWidth', 2);
    end
    hold off; grid off; box off;
    ylabel(['$\bar{z}_s \, | \, \rho = $' num2str(rhoplt(i))]);
    if any(i == [lenrho-1, lenrho])
        xlabel('$s$ (days)');
    end
end
if saveFig
    saveas(gcf, ['allz0Meank_' num2str(lenk) '_' num2str(lenrho)], 'fig');
end

% All worst z0 curves with k for various rho
figure;
for i = 1:lenrho
    subplot(ceil(lenrho/2), 2, i);
    hold on;
    for j = 1:lenk
        plot(z0Worst{i}(j, :), 'Color', cols2{j}, 'LineWidth', 2);
    end
    hold off; grid off; box off;
    ylabel(['$\bar{z}_s \, | \, \rho = $' num2str(rhoplt(i))]);
    if any(i == [lenrho-1, lenrho])
        xlabel('$s$ (days)');
    end
end
if saveFig
    saveas(gcf, ['allz0Worstk_' num2str(lenk) '_' num2str(lenrho)], 'fig');
end

% Pick several k and examine how CI of z_s changes with rho
lex = 4; idex = round(linspace(2, lenk, lex));
% Pick rhoMean will show
lr = 3; idrho = round(linspace(1, lenrho, lr));
% Times to consider - should be fixed across cell
tz = 1:length(z0Mean{1}(1, :)); 

% Figures with colours cols
cols3 = {'r', 'g', 'b'};
figure;
for i = 1:lex
    subplot(ceil(lex/2), 2, i);
    hold on;
    for j = 1:lr
        % CI plots at various rho given k 
        plotCIRaw(tz', z0Mean{idrho(j)}(idex(i), :)', z0Low{idrho(j)}(idex(i), :)',...
            z0High{idrho(j)}(idex(i), :)', cols3{j});
    end
    hold off; grid off; box off;
    ylabel(['$z_s \, | \, k = $' num2str(kplt(idex(i)))]);
    if any(i == [lex-1, lex])
        xlabel('$s$ (days)');
    end
end
if saveFig
    saveas(gcf, ['compCI_' num2str(lenk) '_' num2str(lenrho)], 'fig');
end


% Take specific k and look at rho variations
figure; hold on;
cols4 = {'r', 'g', 'm', 'b'};
for i = 1:lex
    plot(rhoMean, tmeank(idex(i), :), 'Color', cols4{i}, 'LineWidth', 2);
    plot(rhoMean, tworstk(idex(i), :), '--', 'Color', cols4{i}, 'LineWidth', 2);
end
grid off; box off; hold off;
xlabel('$\rho$'); ylabel('$t_{95}$');
if saveFig
    saveas(gcf, ['tdecOverrho_' num2str(lenk) '_' num2str(lenrho)], 'fig');
end




