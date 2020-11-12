% Simple plot from sample trajectory data
clearvars; clc; close all;

% Assumptions and notes
% - summary data from simDieCtrl and extractStats in current folder

% Directory and plotting options
thisDir = cd; [grey1, grey2, cmap] = defaultSet(10);

% Three methods of control/undersampling
ctrlNames = {'uniform', 'subspread', 'superspread'};
ctrlMeth = 1:3; nMeth = length(ctrlMeth);

% Hold all data
data = cell(1, nMeth);
% Limits of all plots (set by trial and error)
x1 = 20; x2 = 250;

% Plot trajectories (averaged across k) for each method
figure;
for i = 1:nMeth
    % Load relevant file (assumed in current dir)
    data{i} = load(join(['statsTraj_' ctrlNames{i}]));
    
    % Time variable should be same across files
    if i == 1
        tday = data{i}.tday;
    end
    
    % Examine mean trajectories and variances with k
    Imeans = cell2mat(data{i}.Imeans'); 
    
    subplot(nMeth, 1, i);
    stairs(tday', Imeans', 'LineWidth', 2);
    grid off; box off;
    xlabel('$s$ (days)', 'FontSize', 18);
    ylabel('E[$I_s$]', 'FontSize', 18);
    title(['Control: ' ctrlNames{i}], 'FontSize', 18);
    xlim([x1 x2]);
end
saveas(gcf, 'exImean', 'fig');

% Plot trajectories (averaged across k) for each method
figure;
for i = 1:nMeth
    % Examine mean trajectories and variances with k
    Ivars = cell2mat(data{i}.Ivars'); 
    
    subplot(nMeth, 1, i);
    stairs(tday', Ivars', 'LineWidth', 2);
    grid off; box off;
    xlabel('$s$ (days)', 'FontSize', 18);
    ylabel('V[$I_s$]', 'FontSize', 18);
    title(['Control: ' ctrlNames{i}], 'FontSize', 18);
    xlim([x1 x2]);
end
saveas(gcf, 'exIvar', 'fig');

% Plot trajectories (averaged across k) for each method
figure;
for i = 1:nMeth
    % Examine mean trajectories and variances with k
    Rmeans = cell2mat(data{i}.Rmeans'); 
    
    subplot(nMeth, 1, i);
    stairs(tday', Rmeans', 'LineWidth', 2);
    grid off; box off;
    xlabel('$s$ (days)', 'FontSize', 18);
    ylabel('E[$R_s$]', 'FontSize', 18);
    title(['Control: ' ctrlNames{i}], 'FontSize', 18);
    xlim([x1 x2]);
end
saveas(gcf, 'exRmean', 'fig');

% Plot trajectories (averaged across k) for each method
figure;
for i = 1:nMeth
    % Examine mean trajectories and variances with k
    Rvars = cell2mat(data{i}.Rvars'); 
    
    subplot(nMeth, 1, i);
    stairs(tday', Rvars', 'LineWidth', 2);
    grid off; box off;
    xlabel('$s$ (days)', 'FontSize', 18);
    ylabel('V[$R_s$]', 'FontSize', 18);
    title(['Control: ' ctrlNames{i}], 'FontSize', 18);
    xlim([x1 x2]);
end
saveas(gcf, 'exRvar', 'fig');