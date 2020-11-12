% Test plot for truncated distributions
clearvars; clc;
close all; tic;

% Assumptions and notes
% - examine various types of truncation with same mean
% - try to match to Lloyd-Smith

% Directory and folder for saving
thisDir = cd; saveFol = 'control type';

% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Main code - compute variance of controls

% Mean R and rho
Rmean = 3; rhoMean = 0.9; 
% Decide number of methods
nMeth = 3; effMean = Rmean*rhoMean;

% Dispersion k of gamma on R
k = logspace(log10(0.5), log10(50), 20); lenk = length(k);  
% Size of R = length(Iday) + zlen
zlen = 100; % no. zeros to append

% No. trajectories to be taken from distribution on Rmean
nTraj = 3000; Rsamp0 = zeros(lenk, nTraj, zlen);
Rsamp1 = Rsamp0; Rsamp2 = Rsamp0;

% Pdf of R for plotting
Rpdf = 0.01:0.01:20; Rprob0 = zeros(lenk, length(Rpdf));
Rprob1 = Rprob0; Rprob2 = Rprob0;

% Draw R samples from gamma distributions based on ctrlMeth
gamMean = zeros(1, lenk); gamVar = gamMean;
% Truncation variables
trMean1 = gamMean; trVar1 = gamMean; trTrunc1 = gamMean;
trMean2 = gamMean; trVar2 = gamMean; trTrunc2 = gamMean;

% Loop through all methods
for j = 1:nMeth
    if j == 1
        % Population wide sampling or control
        for i = 1:lenk
            % Shape-scale gamma
            Rsamp0(i, :, :) = gamrnd(k(i), Rmean*rhoMean/k(i), [nTraj zlen]);
            [gamMean(i), gamVar(i)] = gamstat(k(i), Rmean*rhoMean/k(i));
            % Prob density function of adjusted R
            Rprob0(i, :) = gampdf(Rpdf, k(i), Rmean*rhoMean/k(i));
        end
    elseif j == 2
        % Under-sampling of low reproduction number cases
        for i = 1:lenk
            % Lower-truncated distribution and statistics
            [trDist1, trStat1] = optGamTruncBoth(1, k(i), Rmean*rhoMean/k(i), rhoMean, Rmean, 0.1);
            
            % Store truncation details
            trMean1(i) = trStat1.mean; trVar1(i) = trStat1.var;
            trTrunc1(i) = trDist1.Truncation(1);
            
            % Sample from truncated distribution
            Rsamp1(i, :, :) = random(trDist1, [nTraj zlen]);
            % Prob density function of adjusted R
            Rprob1(i, :) = pdf(trDist1, Rpdf);
        end
    else
        % Control of most infectious - upper truncation of large Rs
        for i = 1:lenk
            % Upper-truncated distribution and statistics
            [trDist2, trStat2] = optGamTruncBoth(2, k(i), Rmean*rhoMean/k(i), rhoMean, Rmean, 10);
            
            % Store truncation details
            trMean2(i) = trStat2.mean; trVar2(i) = trStat2.var;
            trTrunc2(i) = trDist2.Truncation(2);
            
            % Sample from truncated distribution
            Rsamp2(i, :, :) = random(trDist2, [nTraj zlen]);
            % Prob density function of adjusted R
            Rprob2(i, :) = pdf(trDist2, Rpdf);
        end
    end
end

% Compose variance to mean ratios
VM = zeros(nMeth, lenk); VM(1, :) = gamVar./gamMean;
VM(2, :) = trVar1./trMean1; VM(3, :) = trVar2./trMean2;

%% Visualisation and processing

% Make k easier to plot and mean id for saving
k10 = log10(k); meanID = join([num2str(effMean) '_']);
id = strfind(meanID, '.'); meanID(id) = 'x';

% Examine distributions
figure; 
subplot(nMeth, 1, 1);
hold on;
% Middle PDFs
for i = 2:lenk-1
    plot(Rpdf, Rprob0(i, :), 'Color', grey1, 'LineWidth', 2);
end
% Lowest R mean PDF
plot(Rpdf, Rprob0(1, :), 'r', 'LineWidth', 2);
% Laegest R mean PDF
plot(Rpdf, Rprob0(end, :), 'b', 'LineWidth', 2);
grid off; box off; hold off;
xlabel('$R$'); xlim([0 3]);
ylabel('$P(R)$');

subplot(nMeth, 1, 2);
hold on;
% Middle PDFs
for i = 2:lenk-1
    plot(Rpdf, Rprob1(i, :), 'Color', grey1, 'LineWidth', 2);
end
% Lowest R mean PDF
plot(Rpdf, Rprob1(1, :), 'r', 'LineWidth', 2);
% Largest R mean PDF
plot(Rpdf, Rprob1(end, :), 'b', 'LineWidth', 2);
grid off; box off; hold off;
xlabel('$R$'); xlim([0 3]);
ylabel('$P(R)$');

subplot(nMeth, 1, 3);
hold on;
% Middle PDFs
for i = 2:lenk-1
    plot(Rpdf, Rprob2(i, :), 'Color', grey1, 'LineWidth', 2);
end
% Lowest R mean PDF
plot(Rpdf, Rprob2(1, :), 'r', 'LineWidth', 2);
% Largest R mean PDF
plot(Rpdf, Rprob2(end, :), 'b', 'LineWidth', 2);
grid off; box off; hold off;
xlabel('$R$'); xlim([0 3]);
ylabel('$P(R)$');

% Examine means
figure;
subplot(1, 2, 1);
semilogx(k, gamMean, '.-', 'LineWidth', 2, 'MarkerSize', 40);
hold on;
semilogx(k, trMean1, '.-', 'LineWidth', 2, 'MarkerSize', 40);
semilogx(k, trMean2, '.-', 'LineWidth', 2, 'MarkerSize', 40);
grid off; box off; hold off;
xlabel('$k$', 'FontSize', 18);
ylabel('E[$R$]', 'FontSize', 18);
ylim((rhoMean*Rmean)*[0.9 1.1]);

% Examine variances
subplot(1, 2, 2);
semilogx(k, gamVar, '.-', 'LineWidth', 2, 'MarkerSize', 40);
hold on;
semilogx(k, trVar1, '.-', 'LineWidth', 2, 'MarkerSize', 40);
semilogx(k, trVar2, '.-', 'LineWidth', 2, 'MarkerSize', 40);
grid off; box off; hold off;
legend('uniform', 'sub-spread', 'super-spread', 'Location', 'best')
xlabel('$k$', 'FontSize', 18); 
ylabel('V[$R$]', 'FontSize', 18);

% Variance to mean ratio
figure;
plot(k10, VM, 'LineWidth', 2);
h = gca; hL = h.XTickLabel; lenh = length(hL);
for i = 1:lenh
    h.XTickLabel{i} = round(10^str2double(hL{i}), 2, 'significant');
end
xlabel('$k$', 'FontSize', 18);
ylabel('VM[$R$]', 'FontSize', 18);
grid off; box off;
legend('uniform', 'sub-spread', 'super-spread', 'Location', 'best')
cd(saveFol);
saveas(gcf, ['VM_' meanID num2str(nMeth) '_' num2str(lenk)], 'fig');
cd(thisDir);

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% All variable names except R samples
vars = whos; nvar = length(vars);
varNam = cell(1, 1); j = 1;
for i = 1:nvar
   if ~(strncmp(vars(i).name, 'Rsamp', 5))
       % Store this variable
       varNam{j} = vars(i).name;
       j = j + 1;
   end
end

% Save key processed variables
cd(saveFol);
save(join(['ctrlVar_' meanID num2str(nMeth) '_' num2str(lenk)]), varNam{:});
cd(thisDir);
