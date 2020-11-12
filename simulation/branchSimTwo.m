% Simulate Galton-Watson epidemic process elimination
clearvars; clc;
close all; tic;

% Assumptions and notes
% - adds an extra dependence on a generation further
% - uses standard branching process with changing mean
% - declarations are now in absolute time as single 0 ends it

% Folders for saving data
thisDir = cd; saveFol = 'branch sim';
% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Main simulation code

% Num replicate epidemics and control method
M = 5000; ctrlMeth = 1;
% Dispersion k of gamma on R
k = logspace(log10(0.5), 2, 12); lenk = length(k); 
k = 0.5; lenk = 1;

% Switch times and mean reproduction numbers
Rs = [2.5 0.6]; ts = 5; mu = 0.95;
% Define number of days to simulate
nday = 80; tday = 1:nday;

% Name for saving
namstr = [num2str(M) '_' num2str(nday) '_' num2str(lenk) '_' num2str(ctrlMeth)];

%% Simulate branching processes

% Trajectory stats and time of epidemic end
Imeans = cell(1, lenk); Ivars = Imeans; t0s = Imeans; tdecs = Imeans;
z0means = Imeans; z0vars = Imeans; z0mins = Imeans;

% Loop over possible dispersions
for ii = 1:lenk
    % Store variables from each epidemic
    Iday = cell(1, M); t0 = zeros(1, M); z0 = Iday; tdec = t0;
    
    % Sample R by phase of epidemic and control method
    switch(ctrlMeth)
        case 1
            % Population wide sampling or control
            trDist1 = makedist('Gamma', 'a', k(ii), 'b', Rs(1)/k(ii));
            trDist2 = makedist('Gamma', 'a', k(ii), 'b', Rs(2)/k(ii));
        case 2
            % Under-sampling of low reproduction number cases
            [trDist1, ~] = optGamTruncBoth(1, k(ii), Rs(1)/k(ii), 0.5, 2*Rs(1), 0.1);
            [trDist2, ~] = optGamTruncBoth(1, k(ii), Rs(2)/k(ii), 0.5, 2*Rs(2), 0.1);
        case 3
            % Under-sampling of high reproduction number cases
            [trDist1, ~] = optGamTruncBoth(2, k(ii), Rs(1)/k(ii), 0.5, 2*Rs(1), 50);
            [trDist2, ~] = optGamTruncBoth(2, k(ii), Rs(2)/k(ii), 0.5, 2*Rs(2), 10);
    end
    
    for j = 1:M
        % Simulate epidemic scenarios
        Iwarn = 1;
        while Iwarn
            [Iday{j}, Iwarn] = branchSimDieCtrlTwo(nday, ts, trDist1, trDist2);
        end
        % True extinction time
        t0(j) = find(Iday{j}, 1, 'last') + 1;
        % Probability of elimination and declaration time
        z0{j} = exp(-Iday{j}); tdec(j) = find(z0{j} >= mu, 1, 'first');
        disp(['At ' num2str(j) ' of ' num2str(M)]);
    end

    % Convert to matrix and get statistics over time
    Iday = cell2mat(Iday'); Imean = mean(Iday); Ivar = var(Iday);
    z0 = cell2mat(z0'); z0mean = mean(z0); z0var = var(z0); z0min = min(z0);
    
    % Completed M runs at a given k
    disp(['Completed: ' num2str(ii) ' of ' num2str(lenk)]);

    % Statistics over runs
    Imeans{ii} = Imean; Ivars{ii} = Ivar; t0s{ii} = t0; tdecs{ii} = tdec;
    z0means{ii} = z0mean; z0vars{ii} = z0var; z0mins{ii} = z0min;
end


%% Visualisation and plots

% Examine mean trajectories and variances with k
Imeans = cell2mat(Imeans'); Ivars = cell2mat(Ivars');
% Examine statistics of zero cases (end time)
t0mean = cellfun(@mean, t0s); t0var = cellfun(@var, t0s);
t0max = cellfun(@max, t0s); t0s = cell2mat(t0s');

% Strings of k for legend
kstr = cell(1, lenk);
for i = 1:lenk
    kstr{i} = num2str(round(k(i), 2, 'significant'));
end

% Plot the trajectory stats for R and I
figure;
plot(tday, Imeans, 'LineWidth', 2);
grid off; box off;
xlabel('$s$ (days)', 'FontSize', 18); 
ylabel('E[$I_s$]', 'FontSize', 18);
legend(kstr, 'Location', 'best');
cd(saveFol);
saveas(gcf, ['meanITwo_' namstr], 'fig');
cd(thisDir);

figure;
plot(tday, Ivars, 'LineWidth', 2);
grid off; box off;
xlabel('$s$ (days)', 'FontSize', 18); 
ylabel('V[$I_s$]', 'FontSize', 18);
legend(kstr, 'Location', 'best');
cd(saveFol);
saveas(gcf, ['varITwo_' namstr], 'fig');
cd(thisDir);

figure;
semilogx(k, t0mean, 'LineWidth', 2);
hold on;
semilogx(k, t0var, 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$s$ (days)', 'FontSize', 18); 
ylabel('E[$t_0$], V[$t_0$]', 'FontSize', 18);
cd(saveFol);
saveas(gcf, ['t0Two_' namstr], 'fig');
cd(thisDir);

% Data save
cd(saveFol);
save(join(['dataTwo' namstr]));
cd(thisDir);
