% Compute truncation probabilities and VM ratios
function [VM, Rprob0, Rprob1, Rprob2, savStr] = getTruncCtrl(nMeth, k, lenk, Rpdf, effMean, rhoMean, Rmean, saveFol, thisDir, grey1)

% Assumptions and notes
% - works as tesPltTrunc but for batch use
% - modified to compute Rprob as diff of Rcdf

% Probabilities from truncated distributions
Rprob0 = zeros(lenk, length(Rpdf));
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
            [gamMean(i), gamVar(i)] = gamstat(k(i), effMean/k(i));
            % Prob density function of adjusted R
            Rcdf = gamcdf([0 Rpdf], k(i), effMean/k(i));
            % Probabilities as differences of the CDF
            Rprob0(i, :) = diff(Rcdf);
        end
        % Check normalisation
        Rnorm = sum(Rprob0, 2)'; disp(['Normalisation = ' num2str(Rnorm)]);
    elseif j == 2
        % Under-sampling of low reproduction number cases
        for i = 1:lenk
            % Lower-truncated distribution and statistics
            [trDist1, trStat1] = optGamTruncBoth(1, k(i), effMean/k(i), rhoMean, Rmean, 0.1);
            
            % Store truncation details
            trMean1(i) = trStat1.mean; trVar1(i) = trStat1.var;
            trTrunc1(i) = trDist1.Truncation(1);
            
            % Prob density function of adjusted R
            Rcdf = cdf(trDist1, [0 Rpdf]);
            Rprob1(i, :) = diff(Rcdf);
        end
        % Check normalisation
        Rnorm = sum(Rprob1, 2)'; disp(['Normalisation = ' num2str(Rnorm)]);
    else
        % Control of most infectious - upper truncation of large Rs
        for i = 1:lenk
            % Upper-truncated distribution and statistics
            [trDist2, trStat2] = optGamTruncBoth(2, k(i),effMean/k(i), rhoMean, Rmean, 10);
            
            % Store truncation details
            trMean2(i) = trStat2.mean; trVar2(i) = trStat2.var;
            trTrunc2(i) = trDist2.Truncation(2);
            
            % Prob density function of adjusted R
            Rcdf = cdf(trDist2, [0 Rpdf]);
            Rprob2(i, :) = diff(Rcdf);
        end
        % Check normalisation
        Rnorm = sum(Rprob2, 2)'; disp(['Normalisation = ' num2str(Rnorm)]);
    end
end

% Compose variance to mean ratios
VM = zeros(nMeth, lenk); VM(1, :) = gamVar./gamMean;
VM(2, :) = trVar1./trMean1; VM(3, :) = trVar2./trMean2;

% String for saving data/figs
meanID = join([num2str(effMean) '_']);
id = strfind(meanID, '.'); meanID(id) = 'x';
savStr = join([meanID num2str(nMeth) '_' num2str(lenk) '_' num2str(Rmean)]);
% Easier to plot log k
k10 = log10(k); 

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
saveas(gcf, ['VM_' savStr], 'fig');
cd(thisDir);

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
% cd(saveFol);
% saveas(gcf, ['Rdist_' savStr], 'fig');
% cd(thisDir);

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
cd(saveFol);
saveas(gcf, ['meanVar_' savStr], 'fig');
cd(thisDir);
close all;

% Save key processed variables
cd(saveFol);
save(join(['ctrlVar_' savStr]));
cd(thisDir);