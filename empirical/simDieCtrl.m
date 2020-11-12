% Simulate M epidemics at various k but under controls
clearvars; clc;
close all; tic;

% Assumptions and notes
% - uses a single Rmean but epidemic of different k simulated
% - compose elimination probabilities and look at mean and min of z
% - consider theory on E[z] and var(z) against k
% - controlled epidemics now considered

% Folders for saving data
thisDir = cd; saveFol = 'trunc data';
% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Setup M dying epidemics with various k

% Num replicate epidemics and control method
M = 500; ctrlMeth = 2;
% Dispersion k of gamma on R
k = logspace(-1, 2, 12); lenk = length(k); 

% Define number of days to simulate
nday0 = 501; tday0 = 1:nday0;

% Choose scenario and generation time distribution
scenNo = 1; distNo = 2;

% Define all dying epidemic scenarios
scenNam = {'control', 'recovery', 'cascade', 'boom-bust'};
scenChoice = scenNam{scenNo}; disp(['True R scenario: ' scenNam{scenNo}]);

% Define all SI/generation time distributions
distNam = {'exponential', 'gamma', 'delta', 'bimodal'};
distChoice = distNam{distNo}; disp(['True SI scenario: ' distNam{distNo}]);

% Set serial interval parameters
distvals.type = distNo; distvals.omega = 15.3;

% Define SI distribution and R scenario
[distvals.pm, Rs, ts] = setupSIandR(distNo, scenNo);

% Name for saving
namstr = [num2str(M) '_' num2str(scenNo) '_' num2str(distNo) '_'];

%% Simulate epidemic trajectories over k and R

% Variables to interrogate tajectories
Imeans = cell(1, lenk); Ivars = Imeans; 
nday = zeros(1, lenk); Rmeans = Imeans; Rvars = Imeans;

% Loop over possible dispersions
for ii = 1:lenk
    % Store variables from each epidemic
    Iday = cell(1, M); Lam = Iday; Rsamp = Iday; tday = Iday;

    for j = 1:M
        % Simulate epidemic scenarios and truncate
        Iwarn = 1; % ensure no warnings
        while Iwarn
            [Iday{j}, Lam{j}, Rsamp{j}, tday{j}, Iwarn] = epiSimDieCtrl(nday0,...
                ts, Rs, scenNo, distvals, k(ii), ctrlMeth);
        end
        disp(['Run ' num2str(j)]);
    end
    
    % Times should be identical
    tday = cell2mat(tday'); 
    if all(tday(1,:) == tday(end, :))
        tday = unique(tday); nday(ii) = length(tday);
    else
        error('Time stamps inconsistent');
    end

    % Convert to matrix and get statistics over time
    Iday = cell2mat(Iday'); Rsamp = cell2mat(Rsamp'); Lam = cell2mat(Lam');
    Imean = mean(Iday); Ivar = var(Iday);
    Rmean = mean(Rsamp); Rvar = var(Rsamp);
    
    % Completed M runs at a given k
    disp(['Completed: ' num2str(ii) ' of ' num2str(lenk)]);
    % Store this data in a separate file
    cd(saveFol);
    save(join(['sim_' namstr num2str(ii) '.mat']), 'Iday', 'Lam', 'tday', 'Rsamp');
    cd(thisDir);
    
    % Basic testing variables
    Imeans{ii} = Imean; Ivars{ii} = Ivar;
    Rmeans{ii} = Rmean; Rvars{ii} = Rvar;
end

% Key simulation parameters
simPm.M = M; simPm.k = k; simPm.lenk = lenk;
simPm.nday = nday; simPm.distvals = distvals;
simPm.ts = ts; simPm.Rs = Rs; simPm.type = [scenNo distNo];
simPm.Imeans = Imeans; simPm.Ivars = Ivars;

% Also save data with sim parameters
cd(saveFol);
save(join(['setup_' namstr '.mat']), 'simPm');
cd(thisDir);

% Strings of k for legend
kstr = cell(1, lenk);
for i = 1:lenk
    kstr{i} = num2str(round(k(i), 2, 'significant'));
end

% Examine mean trajectories and variances with k
Imeans = cell2mat(Imeans'); Ivars = cell2mat(Ivars');
Rmeans = cell2mat(Rmeans'); Rvars = cell2mat(Rvars');

% Plot the trajectory stats for R and I
figure;
plot(tday, Imeans, 'LineWidth', 2);
grid off; box off;
xlabel('$s$ (days)', 'FontSize', 18); 
ylabel('E[$I_s$]', 'FontSize', 18);
legend(kstr, 'Location', 'best');
cd(saveFol);
saveas(gcf, ['simI_' namstr num2str(lenk)], 'fig');
cd(thisDir);

figure;
plot(tday, Rmeans, 'LineWidth', 2);
grid off; box off;
xlabel('$s$ (days)', 'FontSize', 18);
ylabel('E[$R_s$]', 'FontSize', 18);
legend(kstr, 'Location', 'best');
cd(saveFol);
saveas(gcf, ['simR_' namstr num2str(lenk)], 'fig');
cd(thisDir);

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);




