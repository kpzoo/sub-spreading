% Get true elimination probability given R samples
function [z0Mean, z0High, z0Low, z0Worst, tdecAll] = elimProbCtrl(k, Rmean, mu, rhoMean, ctrlMeth, saveFol, loadFol)

% Assumptions and notes
% - undersampling or control by changing R gamma distribution
% - considers shaping R distribution via ctrlMeth
% - mean of R and rho is fixed and k varies
% - R sampled from a distribution over time and trajectories
% - declarations are relative to tst; when last case was seen
% - examines fixed R vs non-fixed - see projections package

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex', 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Save data and test figs
saveTrue = 1; thisDir = cd;

% Fixed R and no. k to traverse
lenk = length(k);
disp(['Mean of [R, rho] = ' num2str(Rmean) ', ' num2str(rhoMean)]);

%% Load and process empirical data

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

% SI distribution mean and variance
siStat = [15.3, 9.3^2]; % from Djaafara

% Parameters of gamma SI distribution
scalePm = siStat(2)/siStat(1);
shapePm = siStat(1)/scalePm;
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
dom = 0.01:0.01:60;
pdom = gampdf(dom, shapePm, scalePm);
tmaxDec = find(pdom > 10^-4, 1, 'last'); tmaxDec = dom(tmaxDec);

%% Compute renewal elimination probability

% Start checking for epidemic end
tst = nday; nzero = 100;
% Append zeros to Iday (Nishiura used 50)
Iday1 = [Iday, zeros(1, nzero)];

% Times to interrogate for z
tindex = tst+1:length(Iday1);
lenind = length(tindex); tz = 1:lenind;

% No. trajectories to be taken from distribution on Rmean
nTraj = 3000; Rsamp = zeros(lenk, nTraj, zlen);
% Pdf of R for plotting
Rpdf = 0:0.001:10; Rprob = zeros(lenk, length(Rpdf));

% Draw R samples from gamma distributions based on ctrlMeth
gamMean = zeros(1, lenk); gamVar = gamMean;
switch(ctrlMeth)
    case 1
        % Population wide sampling or control
        for i = 1:lenk
            % Shape-scale gamma
            Rsamp(i, :, :) = gamrnd(k(i), Rmean*rhoMean/k(i), [nTraj zlen]);
            [gamMean(i), gamVar(i)] = gamstat(k(i), Rmean*rhoMean/k(i));
            % Prob density function of adjusted R
            Rprob(i, :) = gampdf(Rpdf, k(i), Rmean*rhoMean/k(i));
        end
end

% Elimination probabilities conditioned on R
z0 = zeros(nTraj, lenind); tdec0 = zeros(1, nTraj);
% Capture across R means
z0s = cell(1, lenk); tdec0s = z0s;

% Pre-compute gamma distributions - max length expected
maxlen = length(Iday1)+zlen;
PomegaMax = gampdf(1:maxlen, distvals(1), distvals(2));

% Mean controlling distribution of R
for ii = 1:lenk
    % For t > tst get elimination probability
    for j = 1:nTraj
        for i = 1:lenind
            % True data considered
            Itemp = Iday1(1:tindex(i));
            
            % True elimination probabilities given R
            [z0(j, i), ~] = getProbEliminTrueVecFast(Itemp, squeeze(Rsamp(ii, j, :))', PomegaMax, zlen);
        end
        % Declaration time for a given R
        tdec0(j) = find(z0(j, :) >= mu, 1, 'first');
    end
    % Store results for a given Rmean
    z0s{ii} = z0; tdec0s{ii} = tdec0;
    clear z0 tdec0;
    disp(['Completed ' num2str(ii) ' of ' num2str(lenk)]);
end

% Lower and upper 95% quantile functions
q1 = @(X) quantile(X, 0.025); q2 = @(X) quantile(X, 1-0.025);

% Mean and quantiles of z0s and tdec0s
tdecM = cellfun(@mean, tdec0s);
tdecL = cellfun(q1, tdec0s);
tdecH = cellfun(q2, tdec0s);
z0M = cellfun(@mean, z0s, 'UniformOutput', false);
z0L = cellfun(q1, z0s, 'UniformOutput', false);
z0H = cellfun(q2, z0s, 'UniformOutput', false);

% Worst case declarations
qW = @(X) quantile(X, 0);
z0W = cellfun(qW, z0s, 'UniformOutput', false);
tdecW2 = zeros(1, lenk); tdecM2 = tdecW2;
% Declarations from mean and worst z0 curves
for i = 1:lenk
    ttemp = find(z0W{i} >= mu); ttemp2 = find(z0M{i} >= mu);
    tdecW2(i) = ttemp(1); tdecM2(i) = ttemp2(1);
end
% Other calculation of worst from tdec samples
tdecW = cellfun(@max, tdec0s);
if all(tdecW == tdecW2)
    disp('Worst tdec match');
    idWorst = 1;
end

%% Visualisation and processing

% Moniker for fig saving if k < 1
if Rmean >= 1
    figStr = [num2str(round(Rmean)) '_' num2str(round(10*rhoMean))];
    figStr = join(['Rfixed' figStr]);
else
    figStr = [num2str(round(10*Rmean)) '_' num2str(round(10*rhoMean))];
    figStr = join(['Rfixedpt' figStr]);
end

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
if saveTrue
    cd(saveFol);
    saveas(gcf, ['SIandI_' figStr], 'fig');
    cd(thisDir);
end

% All distributions over R
figure;
hold on;
% Middle PDFs
for i = 2:lenk-1
    plot(Rpdf, Rprob(i, :), 'Color', grey1, 'LineWidth', 2);
end
% Lowest R mean PDF
plot(Rpdf, Rprob(1, :), 'r', 'LineWidth', 2);
% Lareget R mean PDF
plot(Rpdf, Rprob(end, :), 'b', 'LineWidth', 2);
grid off; box off;
xlabel('$R$'); xlim([0 3]);
ylabel('$P(R)$');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['Rdistr_' figStr], 'fig');
    cd(thisDir);
end

% Worst case elimination probabilities and declarations
figure;
subplot(2, 1, 1);
hold on;
% Middle R curves
for i = 2:lenk-1
    stairs(tz, z0W{i}, 'Color', grey1, 'LineWidth', 2);
end
% Lowest R curve
stairs(tz, z0W{1}, 'Color', 'r', 'LineWidth', 2);
% Largest R curve
stairs(tz, z0W{end}, 'Color', 'b', 'LineWidth', 2);
% Maximum from SI
stairs([tmaxDec tmaxDec], [0 1], 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$s$ (days)');
ylabel('$\min z_s$');
xlim([tz(1), round(1.2*tmaxDec)]);
subplot(2, 1, 2);
h = stem(k, tdecW2, 'filled', 'LineWidth', 2);
h.MarkerSize = 10; h.MarkerEdgeColor = 'c'; h.Color = 'c';
hold on;
plot(k, tmaxDec*ones(size(k)), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$k$');
ylabel('$\max t_{95}$');
xlim([k(1) k(end)]); 
h = gca; h.XScale = 'log';
if saveTrue
    cd(saveFol);
    saveas(gcf, ['z0tdecWorst_' figStr], 'fig');
    cd(thisDir);
end

% Mean elimination probabilities and declarations
figure;
subplot(2, 1, 1);
hold on;
% Middle R curves
for i = 2:lenk-1
    stairs(tz, z0M{i}, 'Color', grey1, 'LineWidth', 2);
end
% Lowest R curve
stairs(tz, z0M{1}, 'Color', 'r', 'LineWidth', 2);
% Largest R curve
stairs(tz, z0M{end}, 'Color', 'b', 'LineWidth', 2);
% Maximum from SI
stairs([tmaxDec tmaxDec], [0 1], 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$s$ (days)');
ylabel('$\bar{z}_s$');
xlim([tz(1), round(1.2*tmaxDec)]);
subplot(2, 1, 2);
h = stem(k, tdecM, 'filled', 'LineWidth', 2);
h.MarkerSize = 10; h.MarkerEdgeColor = 'c'; h.Color = 'c';
hold on;
plot(k, tmaxDec*ones(size(k)), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$k$');
ylabel('$\bar{t}_{95}$');
xlim([k(1) k(end)]); ylim([0 tmaxDec+1]);
h = gca; h.XScale = 'log';
if saveTrue
    cd(saveFol);
    saveas(gcf, ['z0tdecMean_' figStr], 'fig');
    cd(thisDir);
end

% Stem plots of all tdecs
figure;
h = stem(k, tdecM, 'filled', 'LineWidth', 2);
h.MarkerSize = 10; h.MarkerEdgeColor = 'b'; h.Color = 'b';
hold on;
h = stem(k, tdecL, 'filled', 'LineWidth', 2);
h.MarkerSize = 10; h.MarkerEdgeColor = grey1; h.Color = grey1;
h = stem(k, tdecH, 'filled', 'LineWidth', 2);
h.MarkerSize = 10; h.MarkerEdgeColor = grey1; h.Color = grey1;
plot(k, tmaxDec*ones(size(k)), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$k$');
ylabel('$t_{95}$');
xlim([k(1) k(end)]); 
h = gca; h.XScale = 'log';
if saveTrue
    cd(saveFol);
    saveas(gcf, ['tdecVary_' figStr], 'fig');
    cd(thisDir);
end

% Consider range of 33% quantiles
id = round(quantile(1:lenk, [0:0.333:1]));
kpltid = k(id); kpltid = round(kpltid, 2, 'significant');

% Capture uncertainty R distribution leaves on z0
figure;
ax = zeros(1, 4);
for i = 1:4
    ax(i) = subplot(2, 2, i);
    plot(tz, z0M{id(i)}, 'Color', 'c', 'LineWidth', 2);
    hold on;
    plotCIRaw(tz', z0M{id(i)}', z0L{id(i)}', z0H{id(i)}', 'c');
    plot(tz, 0.95*ones(size(tz)), 'Color', grey1, 'LineWidth', 2);
    plot([tmaxDec tmaxDec], [0 1.1], 'Color', grey1, 'LineWidth', 2);
    grid off; box off; hold off;
    xlabel('$\Delta s$ (days)');
    ylabel(['$z_s \, | \, k = $' num2str(kpltid(i))]);
    xlim([tz(1), round(1.2*tmaxDec)]);
    ylim([0 1.1]);
end
linkaxes(ax, 'xy');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['CIz_' figStr], 'fig');
    cd(thisDir);
end

% Capture uncertainty R distribution leaves on tdec0
figure;
ax = zeros(1, 4);
for i = 1:4
    ax(i) = subplot(2, 2, i);
    h = histogram(tdec0s{id(i)}, 'Normalization', 'probability');
    hold on;
    plot([tmaxDec tmaxDec], [0 1.01], 'Color', grey1, 'LineWidth', 2);
    grid off; box off; hold off;
    xlabel('$t_{95}$ (days)');
    ylabel('$P(t_{95})$');
    ylim([0 1.01]);
end
linkaxes(ax, 'xy');
if saveTrue
    cd(saveFol);
    saveas(gcf, ['CItdec_' figStr], 'fig');
    cd(thisDir);
end

% Current time of simulation
simDate = string(datetime('now'));

% Outputs of function
z0Mean = cell2mat(z0M'); z0Low = cell2mat(z0L');
z0High = cell2mat(z0H'); z0Worst = cell2mat(z0W');
tdecAll = [tdecM' tdecL' tdecH' tdecW'];

if saveTrue
    cd(saveFol);
    clear Rsamp z0s
    save(join(['elimVaryR_' simDate '_' num2str(lenk) '.mat']));
    cd(thisDir);
end