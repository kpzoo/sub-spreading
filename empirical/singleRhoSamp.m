% Batch across elimination probabilities
clearvars; clc;
close all; tic;

% Assumptions and notes
% - uses a single mean rho but tries different truncations
% - reduced R is what expect to affect future incidence
% - treats control and undersampling in same sense
% - fix Rmean to investigate effect of control via rhoMean
% - either pop wide or indiv based with upper or lower truncation
% - R sampled from a distribution over time and trajectories
% - declarations are relative to tst; when last case was seen
% - incidence curve from R (Ebola) forms the data

% Aditional plotting/partition package
addpath(genpath('/Users/kp10/Documents/MATLAB'));
%addpath(genpath('/Users/kris/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', 18);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Save figures
saveFig = 1; thisDir = cd;
% Folders for loading/saving
loadFol = 'incidence'; saveFol = 'ctrl_single';

%% Main inputs to elimVaryRVaryk and results

% Confidence level for declaration
mu = 0.95; disp(['Confidence = ' num2str(mu)]);
% Mean R and rho
Rmean = 0.9; rhoMean = 0.5; 

% Three methods of control/undersampling
ctrlNames = {'constant', 'lower', 'upper'};
ctrlMeth = [1 3]; lenMeth = length(ctrlMeth);
% Dispersion k of gamma on R
%k = logspace(-2, 1, 10); lenk = length(k);  
k = logspace(log10(0.1), log10(50), 10); lenk = length(k);  

% Load key data from simulations in R
cd(loadFol);
% Incidence curve from R
Iday = csvread("inc.csv", 1,1); Iday = Iday';
nday = length(Iday); tday = 1:nday;
% Dates of cases (datetime array)
dates = readtable("dates.csv"); dates = dates.x;
cd(thisDir);

% Size of R = length(Iday) + zlen
zlen = 100; % no. zeros to append

% SI distribution mean and variance (Djaafara)
siStat = [15.3, 9.3^2]; 
% Parameters of gamma SI distribution
scalePm = siStat(2)/siStat(1); shapePm = siStat(1)/scalePm;
distvals = [shapePm, scalePm];

% Complete distribution over times
Pomega = gampdf(tday, shapePm, scalePm);
% Compute total infectiousness
Lam = zeros(size(Iday));
for i = 2:nday
    % Relevant part of serial distribution
    Lam(i) = sum(Iday(i-1:-1:1).*Pomega(1:i-1));
end

% Examine SI distribution - max should upper bound tdec
dom = 0.01:0.01:60; pdom = gampdf(dom, shapePm, scalePm);
tmaxDec = find(pdom > 10^-4, 1, 'last'); tmaxDec = dom(tmaxDec);

% Input data stucture
inp.Iday = Iday; inp.nday = nday; inp.distvals = distvals; 
inp.zlen = zlen; inp.tmaxDec = tmaxDec; inp.dates = dates;

% Output (results) cells
z0Mean = cell(1, 1); z0High = z0Mean; z0Low = z0Mean;
z0Worst = z0Mean; tdecAll = z0Mean;

% Prob of elimination and declaration stats
for j = 1:lenMeth
    [z0Mean{j}, z0High{j}, z0Low{j}, z0Worst{j}, tdecAll{j}] = elimProbCtrlSing(k, Rmean, mu, rhoMean,...
        ctrlMeth(j), saveFol, inp);
    disp(['Control/sampling method = ' ctrlNames{ctrlMeth(j)}]);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Current time of simulation
simDate = string(datetime('now'));

% Save batch data
cd(saveFol);
save(join(['batchCtrl_' num2str(lenMeth) '_' num2str(lenk) '_' simDate '.mat']));
cd(thisDir);

%% Processing results per Rmean

% For plotting
kplt = round(k, 2, 'significant');
Rmeanplt = round(Rmean, 2, 'significant');
rhoplt = round(rhoMean, 2, 'significant');
% Colours for lenrho
cols = cell(1, lenMeth); 
cols(cellfun(@isempty, cols)) = {grey1};
cols{1} = 'r'; cols{end} = 'b';
% Colours for lenk
cols2 = cell(1, lenk);
cols2(cellfun(@isempty, cols2)) = {grey1};
cols2{1} = 'r'; cols2{end} = 'b';

% Collect mean and worst tdec at every k
tmeank = cell(1, lenk); tworstk = tmeank;
for i = 1:lenk
    ttemp = zeros(1, lenMeth); ttemp2 = ttemp;
    for j = 1:lenMeth
        ttemp(j) = tdecAll{j}(i, 1);
        ttemp2(j) = tdecAll{j}(i, end);
    end
    tmeank{i} = ttemp; tworstk{i} = ttemp2;
end
tmeank = cell2mat(tmeank'); tworstk = cell2mat(tworstk');

% Incidence curve and SI distribution
figure;
subplot(2, 1, 1);
stairs(dates, Iday, 'LineWidth', 2);
grid off; box off;
xlabel('$s$ (days)');
ylabel('$I_s$');
subplot(2, 1, 2);
plot(dom, pdom, 'LineWidth', 2);
grid off; box off;
xlabel('$u$ (days)');
ylabel('$w_u$');
if saveFig
    cd(saveFol);
    saveas(gcf, ['SIandI_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end

% Examine trends against k for each method
figure;
for j = 1:lenMeth
    semilogx(k, tmeank(:, j), 'Color', cols{j}, 'LineWidth', 2);
    if j == 1
        hold on;
    end
    semilogx(k, tworstk(:, j), 'Color', cols{j}, 'LineWidth', 2);
end
xlabel('$k$'); ylabel('$t_{95}$');
grid off; hold off; box off;
if saveFig
    cd(saveFol);
    saveas(gcf, ['tdecOverk_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end

% All mean z0 curves with k for given R and each method
figure;
for i = 1:lenk
    subplot(ceil(lenk/2), 2, i);
    hold on;
    for j = 1:lenMeth
        plot(z0Mean{j}(i, :), 'Color', cols{j}, 'LineWidth', 2);
    end
    hold off; grid off; box off;
    ylabel(['$\bar{z}_s \, | \, k = $' num2str(kplt(i))]);
    if any(i == [lenk-1, lenk])
        xlabel('$s$ (days)');
    end
end
if saveFig
    cd(saveFol);
    saveas(gcf, ['allz0MeanRho_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end

% All worst z0 curves with k for given R and each method
figure;
for i = 1:lenk
    subplot(ceil(lenk/2), 2, i);
    hold on;
    for j = 1:lenMeth
        plot(z0Worst{j}(i, :), 'Color', cols{j}, 'LineWidth', 2);
    end
    hold off; grid off; box off;
    ylabel(['$z_{s, min} \, | \, k = $' num2str(kplt(i))]);
    if any(i == [lenk-1, lenk])
        xlabel('$s$ (days)');
    end
end
if saveFig
    cd(saveFol);
    saveas(gcf, ['allz0WorstRho_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end

% All mean z0 curves with k for various rho
figure;
for i = 1:lenMeth
    subplot(ceil(lenMeth/2), 2, i);
    hold on;
    for j = 2:lenk-1
        plot(z0Mean{i}(j, :), 'Color', cols2{j}, 'LineWidth', 2);
    end
    plot(z0Mean{i}(1, :), 'Color', cols2{1}, 'LineWidth', 2);
    plot(z0Mean{i}(end, :), 'Color', cols2{end}, 'LineWidth', 2);
    hold off; grid off; box off;
    ylabel(['$\bar{z}_s \, | \, \rho = $' num2str(rhoplt)]);
    if any(i == [lenMeth-1, lenMeth])
        xlabel('$s$ (days)');
    end
end
if saveFig
    cd(saveFol);
    saveas(gcf, ['allz0Meank_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end

% All worst z0 curves with k for various rho
figure;
for i = 1:lenMeth
    subplot(ceil(lenMeth/2), 2, i);
    hold on;
    for j = 2:lenk-1
        plot(z0Worst{i}(j, :), 'Color', cols2{j}, 'LineWidth', 2);
    end
    plot(z0Worst{i}(1, :), 'Color', cols2{1}, 'LineWidth', 2);
    plot(z0Worst{i}(end, :), 'Color', cols2{end}, 'LineWidth', 2);
    hold off; grid off; box off;
    ylabel(['$\bar{z}_s \, | \, \rho = $' num2str(rhoplt)]);
    if any(i == [lenMeth-1, lenMeth])
        xlabel('$s$ (days)');
    end
end
if saveFig
    cd(saveFol);
    saveas(gcf, ['allz0Worstk_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end

% Pick several k and examine how CI of z_s changes with rho
lex = 4; idex = round(linspace(2, lenk, lex));
% Pick rhoMean will show
lr = 3; idrho = round(linspace(1, lenMeth, lr));
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
    cd(saveFol);
    saveas(gcf, ['compCI_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end

% Look at tdec mean and worst for each method
figure;
subplot(2, 1, 1); 
semilogx(k, tmeank, '.-', 'LineWidth', 2, 'MarkerSize', 40);
hold on; semilogx(k, tmaxDec*ones(size(k)), 'k--', 'LineWidth', 2);
grid off; box off; hold off;
xlabel('$k$'); ylabel('$\bar{t}_{95}$');
subplot(2, 1, 2); 
semilogx(k, tworstk, '.-', 'LineWidth', 2, 'MarkerSize', 40);
hold on; semilogx(k, tmaxDec*ones(size(k)), 'k--', 'LineWidth', 2);
grid off; box off; hold off;
xlabel('$k$'); ylabel('$\max \, t_{95}$');
if saveFig
    cd(saveFol);
    saveas(gcf, ['tdecMeth_' num2str(lenk) '_' num2str(lenMeth)], 'fig');
    cd(thisDir);
end




