% Deeper look at truncated distributions
clearvars; clc;
close all; tic;

% Assumptions and notes
% - looks at how truncation control changes low k distributions

% Folders for saving data
thisDir = cd; saveFol = 'branch sim';
% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% View truncated distributions

% Dispersion k and no. samples
k = 0.5; nSamp = 10^6;
% Control degradation and control methods
ctrlMeth = 1:3; actrl = 1;
% Control names
ctrlName = {'uniform', 'sub-spreading', 'super-spreading'};

% Choose case
exHigh = 1; 
if exHigh
    % Large R (mean) case
    Rs = 5; a = 2; b = 18;
else
    % Small R (mean) case
    Rs = 0.5; a = 0.2; b = 1.8;
end
disp(['Rs = ' num2str(Rs)]);

% Effective R
effR = Rs*actrl; disp(['Effective R: ' num2str(effR)])

% Samples from all distributions
xUnif = zeros(1, nSamp); xSub = xUnif; xSup = xUnif; 

% Domain of R for plotting 
%Rdom = 0.0001:0.0001:100; lenR = length(Rdom);
Rdom = linspace(0.0001, 100, 10^5); lenR = length(Rdom);
% CDF variables
Runif = zeros(size(Rdom)); Rsub = Runif; Rsup = Runif; 

% Construct sample distribution for R
for j = ctrlMeth
    switch(j)
        case 1
            % Population wide sampling or control
            trUnif = makedist('Gamma', 'a', k, 'b', Rs/k);
            % Samples and CDF
            xUnif = random(trUnif, [1 nSamp]); Runif = cdf(trUnif, Rdom); 
            % PDF and VM ratios
            Punif = pdf(trUnif, Rdom); VMunif = var(trUnif)/mean(trUnif);
            
        case 2
            % Under-sampling of low reproduction number cases
            [trSub, ~] = optGamTruncTest(1, k, Rs/k, actrl, Rs/actrl, a);
            % Samples and CDF
            xSub = random(trSub, [1 nSamp]); Rsub = cdf(trSub, Rdom); 
            % PDF and VM ratios
            Psub = pdf(trSub, Rdom); VMsub = var(trSub)/mean(trSub); 
            
        case 3
            % Under-sampling of high reproduction number cases
            [trSup, ~] = optGamTruncTest(2, k, Rs/k, actrl, Rs/actrl, b);
            % Samples and CDF
            xSup = random(trSup, [1 nSamp]); Rsup = cdf(trSup, Rdom); 
            % PDF and VM ratios
            Psup = pdf(trSup, Rdom); VMsup = var(trSup)/mean(trSup);
    end
end

% Display ranking of variance to means
VMset = [VMunif VMsub VMsup]; disp(['VM [unif sub sup] = ' num2str(VMset)]);
% Scale parameters of each distribution
bset = [trUnif.b, trSub.b, trSup.b]; disp(['Scale (b): ' num2str(bset)]);
% Mean of all distributions 
mset = [mean(trUnif) mean(trSub) mean(trSup)]; disp(['E[R] [unif sub sup] = ' num2str(mset)]); 

%% Visualisation

% Compare CDFs of distributions
figure;
semilogx(Rdom, Runif, 'Color', 'g', 'LineWidth', 2);
hold on;
semilogx(Rdom, Rsup, 'Color', 'r', 'LineWidth', 2);
semilogx(Rdom, Rsub, 'Color', 'b', 'LineWidth', 2);
semilogx([a a], [0 1], '--', 'Color', grey1, 'LineWidth', 2);
semilogx([b b], [0 1], '--', 'Color', grey1, 'LineWidth', 2);
box off; grid off; hold off;
xlim([0.01 Rdom(end)]);
xlabel('$R$', 'FontSize', 18);
ylabel('F($R$)', 'FontSize', 18);
title(['$\rho \mu$ = ' num2str(Rs)], 'FontSize', 18);

% Bar plots of key stats
figure;
subplot(2, 2, 1);
% Means 
h = bar(mset); h.FaceColor = 'b'; h.FaceAlpha = 0.2;
box off; grid off; 
xlabel('control type', 'FontSize', 18);
ylabel('E[$R$]', 'FontSize', 18);
subplot(2, 2, 2);
% VM ratios 
h = bar(VMset); h.FaceColor = 'b'; h.FaceAlpha = 0.2;
box off; grid off; 
xlabel('control type', 'FontSize', 18);
ylabel('VM[$R$]', 'FontSize', 18);
subplot(2, 2, 3);
% Scale parameters
h = bar(bset); h.FaceColor = 'b'; h.FaceAlpha = 0.2;
box off; grid off; 
xlabel('control type', 'FontSize', 18);
ylabel('scale', 'FontSize', 18);
h = subplot(2, 2, 4);
% Parameter settings
%uitable('Data', [Rs a b k]', 'RowName', {'R', 'a', 'b', 'k'}, 'Units', 'Normalized', 'Position', h.Position);
text(0.1, 0.1, ['Mean $R$ = ' num2str(Rs)], 'FontSize', 18); 
text(0.1, 0.3, ['Lower limit $a$ = ' num2str(a)], 'FontSize', 18); 
text(0.1, 0.5, ['Upper limit $b$ = ' num2str(b)], 'FontSize', 18); 
text(0.1, 0.7, ['Dispersion $k$ = ' num2str(k)], 'FontSize', 18); 
h.XColor = 'none'; h.YColor = 'none';
title('Simulation parameters', 'FontSize', 18);

% Corresponding PDFs on low scale
figure;
subplot(1, 3, 1);
plot(Rdom(Rdom <= a), Punif(Rdom <= a), 'Color', 'g', 'LineWidth', 2);
hold on;
loglog(Rdom(Rdom <= a), Psup(Rdom <= a), 'Color', 'r', 'LineWidth', 2);
loglog(Rdom(Rdom <= a), Psub(Rdom <= a), 'Color', 'b', 'LineWidth', 2);
box off; grid off; hold off;
xlim([0.01 a]); 
xlabel('$R$', 'FontSize', 18);
ylabel('P($R$)', 'FontSize', 18);
title(['$\rho \mu$ = ' num2str(Rs)], 'FontSize', 18);
% Corresponding PDFs on high scale
subplot(1, 3, 2);
plot(Rdom(Rdom >= b), Punif(Rdom >= b), 'Color', 'g', 'LineWidth', 2);
hold on;
loglog(Rdom(Rdom >= b), Psup(Rdom >= b), 'Color', 'r', 'LineWidth', 2);
loglog(Rdom(Rdom >= b), Psub(Rdom >= b), 'Color', 'b', 'LineWidth', 2);
box off; grid off; hold off;
xlim([b Rdom(end)]); 
xlabel('$R$', 'FontSize', 18);
ylabel('P($R$)', 'FontSize', 18);
title(['$\rho \mu$ = ' num2str(Rs)], 'FontSize', 18);
% Corresponding PDFs on mid scale
subplot(1, 3, 3);
plot(Rdom(Rdom >= a & Rdom <= b), Punif(Rdom >= a & Rdom <= b), 'Color', 'g', 'LineWidth', 2);
hold on;
loglog(Rdom(Rdom >= a & Rdom <= b), Psup(Rdom >= a & Rdom <= b), 'Color', 'r', 'LineWidth', 2);
loglog(Rdom(Rdom >= a & Rdom <= b), Psub(Rdom >= a & Rdom <= b), 'Color', 'b', 'LineWidth', 2);
box off; grid off; hold off;
xlim([a b]); 
xlabel('$R$', 'FontSize', 18);
ylabel('P($R$)', 'FontSize', 18);
title(['$\rho \mu$ = ' num2str(Rs)], 'FontSize', 18);

%% Examine scale behaviour of gamma

% Mean reproduction number
Rm = [4 0.5]; lenR = length(Rm);








