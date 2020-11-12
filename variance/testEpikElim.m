% Examine variance/mean z_s for epidemics of various k
clearvars; clc;
close all; tic;

% Assumptions and notes
% - assumes R constant into future
% - uses a single Rmean but epidemic of different k simulated
% - compose elimination probabilities and look at mean and min of z
% - consider theory on E[z] and var(z) against k
% - each type of control applied in turn to each trajectory
% - does not add extra 0s, just run sim for longer

% Directory and confidence
thisDir = cd; mu = 0.95;
% Folders for loading and saving
loadFol = 'sim data'; saveFol = 'sim process';

% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Load data for each k - incidence and R samples

% Num of replicate trajectories
M = 2000; Mstr = num2str(M);
% Update load folder name
loadFol = join([loadFol '/' Mstr]);

% Setup data including range of k
cd(loadFol);
setFile = dir('setup*');
% Simulation setup data
setData = load(setFile.name);
setData = setData.simPm;
cd(thisDir);

% Key simulation variables
k = setData.k; distvals = setData.distvals; 
lenk = setData.lenk; nday = setData.nday(1);

% Constant R assumed into future
Rlast = setData.Rs(end);
% First R and switch time (assumes sinlge switch)
Rfirst = setData.Rs(1); tsw = setData.ts;

% Psuedo zeros to add for computation of zs
zlen = 100; maxlen = zlen + nday;
% Only using gamma SI so convert to parameters
shapePm = distvals.pm; scalePm = distvals.omega/shapePm;
% For extra length define SI distribution
PomegaMax = gampdf(1:maxlen, shapePm, scalePm);

%% Compute zs for every epidemic

% Variance and mean of I(t), R(t) versus k
Rkmean = zeros(nday, lenk); Rkvar = Rkmean;
Ikmean = Rkmean; Ikvar = Rkmean;

% Time of last 0 and declarations
tsts = zeros(lenk, M); tdecs = tsts;
% Elimination probabilities and statistics
zs = cell(lenk, M); 

for ii = 1:lenk
   % Load data for specific k
   cd(loadFol);
   % File name
   fileProp = dir(join(['*_' num2str(ii) '.mat']));
   if length(fileProp) > 1
       error('File naming issues');
   end
   fileProp = fileProp.name;
   % Extract trajectory data
   data = load(fileProp);
   cd(thisDir);
   
   % Every trajectory at this k
   Iday = data.Iday; Rsamp = data.Rsamp; Lam = data.Lam;
   
   % Compute variance/mean of I(t) and R(t) with k
   Ikmean(:, ii) = mean(Iday); Ikvar(:, ii) = var(Iday);
   Rkmean(:, ii) = mean(Rsamp); Rkvar(:, ii) = var(Rsamp);
   
   for j = 1:M
      % Last time of non-zero incidence
      tst = find(Iday(j, :), 1, 'last');
      % Times z computed over
      tindex = tst+1:nday; lenind =  length(tindex);
      
      % Sample of R for the zlen computations
      Rzlen = gamrnd(k(ii), Rlast/k(ii), [1 zlen]);
      
      % Probability of elimination
      z = zeros(1, lenind);
      for i = 1:lenind
          [z(i), ~] = getProbEliminTrueVecFast(Iday(j, 1:tindex(i)), Rzlen, PomegaMax, zlen);
      end
      zs{ii, j} = z;
      
      % Declaration time
      tdec = find(z >= mu, 1, 'first');
      if isempty(tdec)
          disp(['Epidemic did not end, [ii j] = ' num2str(ii) ' ' num2str(j)]);
          tdec = -1;
      end
      % Store times for given trajectory
      tdecs(ii, j) = tdec; tsts(ii, j) = tst;
   end
   disp(['Completed ' num2str(ii) ' of ' num2str(lenk)]);
end

% Statistics of declaration times
tdecM = zeros(1, lenk); tdecV = tdecM; tstM = tdecM;
tstV = tdecM; tdecQ = zeros(3, lenk); tstQ = tdecQ; 
for i = 1:lenk
    tdecM(i) = mean(tdecs(i, tdecs(i, :) ~= -1));
    tdecV(i) = var(tdecs(i, tdecs(i, :) ~= -1));
    tdecQ(:, i) = quantile(tdecs(i, tdecs(i, :) ~= -1), [0.025, 0.5, 0.975]);
    tstQ(:, i) = quantile(tsts(i, tdecs(i, :) ~= -1), [0.025, 0.5, 0.975]);
    tstM(i) = mean(tsts(i, tdecs(i, :) ~= -1));
    tstV(i) = var(tsts(i, tdecs(i, :) ~= -1));
end

% Variance to mean ratios
VMtdec = tdecV./tdecM; k10 = log10(k);
VM_R = Rkvar./Rkmean; VM_I = Ikvar./Ikmean;

% Elimination probabilities at various thresholds
thresh = 0.05:0.05:mu; nthresh = length(thresh);
tthreshk = cell(1, lenk);

for ii = 1:lenk
    % Threshold values for this k
    tthresh = zeros(M, nthresh);
    for j = 1:M
        % Extract z curve
        z = zs{ii, j};
        % Get various threshold times
        if max(z) >= mu
            for i = 1:nthresh
                tthresh(j, i) = find(z >= thresh(i), 1, 'first');
            end
        else
            % Epidemic did not end
            thresh(j, :) = -1;
        end
    end
    % Store all thresholds
    tthreshk{ii} = tthresh;
end

%% Visualisations and processing

% The variance and mean of declaration times
figure;
plot(k10, VMtdec, '.-', 'LineWidth', 2, 'MarkerSize', 40);
ylabel('VM[$t_{95}$]', 'FontSize', 18);
xlabel('$\log_{10} k$', 'FontSize', 18); box off; grid off;

% Compare I and R variance to mean ratios
figure;
subplot(2, 1, 1);
meshz(k10, 1:nday, log10(VM_R));
view([75 25]); grid off; box off;
xlabel('$k$', 'FontSize', 18);
ylabel('$s$ (days)', 'FontSize', 18);
ylim([1 nday]); 
h = gca; hL = h.XTickLabel; lenh = length(hL);
for i = 1:lenh
    h.XTickLabel{i} = 10^str2double(hL{i});
end
hL = h.ZTickLabel; lenh = length(hL);
for i = 1:lenh
    h.ZTickLabel{i} = 10^str2double(hL{i});
end
zlabel('VM[$R_s$]', 'FontSize', 18);

subplot(2, 1, 2);
s = meshz(k10, 1:nday, log10(VM_I));
view([75 25]); grid off; box off;
xlabel('$k$', 'FontSize', 18);
ylabel('$s$ (days)', 'FontSize', 18);
ylim([1 nday]);
h = gca; hL = h.XTickLabel; lenh = length(hL);
for i = 1:lenh
    h.XTickLabel{i} = 10^str2double(hL{i});
end
hL = h.ZTickLabel; lenh = length(hL);
for i = 1:lenh
    h.ZTickLabel{i} = 10^str2double(hL{i});
end
zlabel('VM[$I_s$]', 'FontSize', 18);
cd(saveFol);
saveas(gcf, ['VM_R_I_' Mstr '_' num2str(lenk)], 'fig');
cd(thisDir);

% Quantiles of declaration times given last 0
figure;
plotCI2(k10', tdecM', tdecM'-tdecQ(1,:)', tdecQ(3,:)'-tdecM', 'b', 1);
grid off; box off;
ylabel('$\Delta t_{95}$ (days)', 'FontSize', 18);
xlabel('$\log_{10} k$', 'FontSize', 18);
cd(saveFol);
saveas(gcf, ['tdec_' Mstr '_' num2str(lenk)], 'fig');
cd(thisDir);

% Quantiles of time of last 0
figure;
plotCI2(k10', tstM', tstM'-tstQ(1,:)', tstQ(3,:)'-tstM', 'b', 1);
grid off; box off;
ylabel('$t_{0}$ (days)', 'FontSize', 18);
xlabel('$\log_{10} k$', 'FontSize', 18);
cd(saveFol);
saveas(gcf, ['tst_' Mstr '_' num2str(lenk)], 'fig');
cd(thisDir);

% Mean and variance of R
figure;
subplot(2, 1, 1);
semilogy(1:nday, Rkmean, 'LineWidth', 2);
ylabel('E[$R_s$]', 'FontSize', 18); 
box off; grid off;
subplot(2, 1, 2);
semilogy(1:nday, Rkvar, 'LineWidth', 2);
ylabel('V[$R_s$]', 'FontSize', 18); 
box off; grid off;
cd(saveFol);
saveas(gcf, ['Rstat_' Mstr '_' num2str(lenk)], 'fig');
cd(thisDir);

% Mean and variance of I
figure;
subplot(2, 1, 1);
semilogy(1:nday, Ikmean, 'LineWidth', 2);
ylabel('E[$I_s$]', 'FontSize', 18); 
box off; grid off;
subplot(2, 1, 2);
semilogy(1:nday, Ikvar, 'LineWidth', 2);
ylabel('V[$I_s$]', 'FontSize', 18); 
box off; grid off;
cd(saveFol);
saveas(gcf, ['Istat_' Mstr '_' num2str(lenk)], 'fig');
cd(thisDir);

% 

% Save key processed variables
cd(saveFol);
save(join(['proc_' Mstr '_' num2str(lenk)]), 'tstM', 'tstQ', 'tdecM', 'tdecQ', 'tstQ',...
    'zs', 'Ikmean', 'Ikvar', 'Rkmean', 'Rkvar');
cd(thisDir);

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
