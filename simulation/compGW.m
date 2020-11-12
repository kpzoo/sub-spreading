% Process declaration times from Galton-Watson
clearvars; clc;
close all; tic;

% Assumptions and notes
% - loads simulated Galton-Watson data
% - declarations are now in absolute time as single 0 ends it

% Folders for loading/saving data
thisDir = cd; loadFol = 'branch sim';
% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Load data and compare at smallest k

% Possible control methods
ctrlMeth = 1:3; nMeth = length(ctrlMeth);
% Store main data
mainData = cell(1, nMeth); meth = zeros(1, nMeth);

% Data for each method
cd(loadFol);
% All files of interest
file = dir('GW*'); 
if length(file) ~= nMeth
    error('Wrong number of input files');
end

% For each file get data of interest
for i = 1:nMeth
    % Part of name with method
    filnam = file(i).name;
    id = regexp(filnam, '.mat'); id = id-1;
    meth(i) = str2double(filnam(id));
    
    % Extract data but place indices in order
    mainData{meth(i)} = load(filnam, 'Ivars', 'Imeans', 'tdecs', 'tdecmean', 'tdecvar', 'tdecmax', 'ctrlMeth');
    
    % Setup data
    if meth(i) == 1
        setupData = load(filnam, 'k', 'M', 'lenk', 'nday');
        nday = setupData.nday; M = setupData.M;
        k = setupData.k; lenk = setupData.lenk;
    end
end
cd(thisDir);

% ID of k index to consider
kid = 2; kchoice = k(kid);
disp(['Examining k = ' num2str(kchoice)]);

% Extract declaration times and mean/var of incidence for kid
tdeck = cell(1, nMeth); Imeank = tdeck; Ivark = tdeck;
for i = 1:nMeth
    tdeck{i} = mainData{i}.tdecs{kid};
    Imeank{i} = mainData{i}.Imeans(kid, :);
    Ivark{i} = mainData{i}.Ivars(kid, :);
end

% Get all tdec variance and means across 
tdecvar = cell(1, nMeth); tdecmean = tdecvar;
for i = 1:nMeth
   tdecvar{i} = mainData{i}.tdecvar;
   tdecmean{i} = mainData{i}.tdecmean;
end

%% Visualisations

% Colours
cols{1} = 'g'; cols{2} = 'b'; cols{3} = 'r';

% CDF (empirical) of declaration times and VM ratios
xcdf = cell(1, nMeth); Fcdf = xcdf; VM = xcdf;
for i = 1:nMeth
    % CDF of tdec (95)
    [Fcdf{i}, xcdf{i}] = ecdf(tdeck{i});
    % Variance to mena ratio
    VM{i} = Ivark{i}./Imeank{i};
end
% Plot for each method in turn
figure; hold on;
for i = 1:nMeth
    stairs(xcdf{i}, Fcdf{i}, 'Color', cols{i}, 'LineWidth', 2);
end
hold off; grid off; box off;
xlabel('$t_{95}$ (generations)', 'FontSize', 18); 
ylabel('F($t_{95}$)', 'FontSize', 18);
cd(loadFol);
saveas(gcf, 'cdf_tdec', 'fig');
cd(thisDir);

% Histogram of end times
figure;
% Base case of no truncation
h = histogram(tdeck{1}, 'Normalization', 'probability'); 
hold on; h.FaceColor = cols{1}; 
h.FaceAlpha = 0.2; h.EdgeAlpha = 0.1;
% Other truncations
for i = 2:nMeth
    h = histogram(tdeck{i}, 'Normalization', 'probability');
    h.FaceAlpha = 0.2; h.EdgeAlpha = 0.1;
    h.FaceColor = cols{i}; 
end
hold off; grid off; box off;
xlabel('$t_{95}$ (generations)', 'FontSize', 18); 
ylabel('P($t_{95}$)', 'FontSize', 18);
cd(loadFol);
saveas(gcf, 'hist_tdec', 'fig');
cd(thisDir);

% Mean values of incidence as a check
figure; hold on;
for i = 1:nMeth
    plot(1:nday, Imeank{i}, 'Color', cols{i}, 'LineWidth', 2);
end
hold off; grid off; box off;
xlabel('$s$ (generations)', 'FontSize', 18); 
ylabel('E[$I_{s}$]', 'FontSize', 18);

% VM of incidence 
figure; hold on;
for i = 1:nMeth
    plot(1:nday, VM{i}, 'Color', cols{i}, 'LineWidth', 2);
end
hold off; grid off; box off;
xlabel('$s$ (generations)', 'FontSize', 18); 
ylabel(['VM$[I_{s}] \, | \, k$ = ' num2str(kchoice)], 'FontSize', 18);
cd(loadFol);
saveas(gcf, 'VM_I', 'fig');
cd(thisDir);

% Save the VM data for plotting later
save('VM_GW.mat', 'nday', 'VM', 'cols', 'kchoice');

% Mean and variance of the declaration times
figure;
subplot(2, 1, 1);
for i = 1:nMeth
    semilogx(k, tdecvar{i}, 'Color', cols{i}, 'LineWidth', 2);
    if i == 1
        hold on;
    end
end
hold off; grid off; box off;
xlabel('$k$', 'FontSize', 18); 
ylabel('V[$t_{95}$]', 'FontSize', 18);
xlim([0.1 1]);
subplot(2, 1, 2);
for i = 1:nMeth
    semilogx(k, tdecmean{i}, 'Color', cols{i}, 'LineWidth', 2);
    if i == 1
        hold on;
    end
end
hold off; grid off; box off;
xlim([0.1 1]);
xlabel('$k$', 'FontSize', 18); 
ylabel('E[$t_{95}$]', 'FontSize', 18);
cd(loadFol);
saveas(gcf, 'tdecstats', 'fig');
cd(thisDir);