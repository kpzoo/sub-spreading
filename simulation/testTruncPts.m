% Deeper look at truncation points of control
clearvars; clc;
close all; tic;

% Assumptions and notes
% - look at how variations in a and b matter given mean R

% Folders for saving data
thisDir = cd; saveFol = 'test dist';
% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% View truncated distributions

% Dispersion k 
k = 0.5; 
% Control degradation and control methods
ctrlMeth = 1:3; actrl = 1;
% Control names
ctrlName = {'uniform', 'sub-spreading', 'super-spreading'};

% Choose means and a (lower), b (higher) truncation
exHigh = 0; 
if exHigh
    % Large R (mean) case
    Rs = 5; lena = 50;
    a = linspace(0.01, 2, lena);
    b = linspace(15.5, 80, lena);
else
    % Small R (mean) case
    Rs = 0.5; lena = 50;
    a = linspace(0.01, 0.2, lena);
    b = linspace(3, 10, lena);
end
disp(['Rs = ' num2str(Rs)]);

% Name for saving
namstr = [num2str(actrl) '_' num2str(exHigh)];

% Effective R 
effR = Rs*actrl; disp(['Effective R: ' num2str(effR)])
% Domain of R for plotting 
Rdom = 0.0001:0.0001:100; lenR = length(Rdom);

% CDF variables
Runif = cell(1, lena); Rsub = Runif; Rsup = Runif; 
% VM variables
VMunif = zeros(1, lena); VMsub = VMunif; VMsup = VMunif; 
% Altered scale parameters
sc_unif = zeros(1, lena); sc_sub = sc_unif; sc_sup = sc_unif;
munif = zeros(1, lena); msub = munif; msup = munif;

% Construct sample distribution for R
for i = 1:lena
    for j = ctrlMeth
        switch(j)
            case 1
                % Population wide sampling or control
                trUnif = makedist('Gamma', 'a', k, 'b', Rs/k);
                % CDF and VM ratios
                Runif{i} = cdf(trUnif, Rdom); VMunif(i) = var(trUnif)/mean(trUnif);
                % New scale parameters and means
                sc_unif(i) = trUnif.b; munif(i) = mean(trUnif);
            case 2
                % Under-sampling of low reproduction number cases
                [trSub, ~] = optGamTruncTest(1, k, Rs/k, actrl, Rs/actrl, a(i));
                % CDF and VM ratios
                Rsub{i} = cdf(trSub, Rdom); VMsub(i) = var(trSub)/mean(trSub);
                % New scale parameters and means
                sc_sub(i) = trSub.b; msub(i) = mean(trSub);
                
            case 3
                % Under-sampling of high reproduction number cases
                [trSup, ~] = optGamTruncTest(2, k, Rs/k, actrl, Rs/actrl, b(i));
                % CDF and VM ratios
                Rsup{i} = cdf(trSup, Rdom); VMsup(i) = var(trSup)/mean(trSup);
                % New scale parameters and means
                sc_sup(i) = trSup.b; msup(i) = mean(trSup);
        end
    end
end

%% Visualisation

% VMs against truncation points
figure;
subplot(2, 1, 1);
h = bar(b, VMsup);
h.FaceColor = 'r'; h.FaceAlpha = 0.2; 
hold on;
plot(b, VMunif, 'Color', 'g', 'LineWidth', 2);
box off; grid off; hold off;
xlabel('upper limit $b$', 'FontSize', 18);
ylabel(['VM[$R$]$\, | \,\mu = $' num2str(Rs)], 'FontSize', 18);
subplot(2, 1, 2);
h = bar(a, VMsub);
h.FaceColor = 'b'; h.FaceAlpha = 0.2; 
hold on;
plot(a, VMunif, 'Color', 'g', 'LineWidth', 2);
box off; grid off; hold off;
xlabel('lower limit $a$', 'FontSize', 18);
ylabel(['VM[$R$]$\, | \,\mu = $' num2str(Rs)], 'FontSize', 18);
cd(saveFol);
saveas(gcf, ['VMR_' namstr], 'fig');
cd(thisDir);

% VMs against truncation points
figure;
subplot(2, 1, 1);
h = bar(b, sc_sup);
h.FaceColor = 'r'; h.FaceAlpha = 0.2; 
hold on;
plot(b, sc_unif, 'Color', 'g', 'LineWidth', 2);
box off; grid off; hold off;
xlabel('upper limit $b$', 'FontSize', 18);
ylabel(['Scaled $\theta | \,\mu = $' num2str(Rs)], 'FontSize', 18);
subplot(2, 1, 2);
h = bar(a, sc_sub);
h.FaceColor = 'b'; h.FaceAlpha = 0.2; 
hold on;
plot(a, sc_unif, 'Color', 'g', 'LineWidth', 2);
box off; grid off; hold off;
xlabel('lower limit $a$', 'FontSize', 18);
ylabel(['Scaled $\theta | \,\mu = $' num2str(Rs)], 'FontSize', 18);
cd(saveFol);
saveas(gcf, ['Scale_' namstr], 'fig');
cd(thisDir);


% Data save
cd(saveFol);
clear('Rsub', 'Rsup', 'Runif', 'Rdom');
save(join(['exHigh' namstr]));
cd(thisDir);








