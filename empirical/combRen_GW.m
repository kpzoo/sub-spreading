% Combine plots of VM from GW and renewal model
clearvars; clc; close all;

% Assumptions and notes
% - all data in ctrlsim folder, place GW data there

% Folders for saving data
thisDir = cd; loadFol = 'ctrlsim';
% Directory of some main code and plotting options
cd ..; mainDir = cd; mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);


% Three methods of control/undersampling
ctrlNames = {'uniform', 'subspread', 'superspread'};
ctrlMeth = 1:3; nMeth = length(ctrlMeth);

% Extract all VM data from renewal sims
cd(loadFol);
% All renewal files
files = dir('sim_*'); nFil = length(files);
% Should match control methods
if nFil ~= nMeth
    error('File count seems incorrect');
end

% VM storage variables
VM_Is = cell(1, nMeth); 
for i = 1:nFil
   % Extract VM data 
   data = load(files(i).name, 'VM_I', 'tday');
   VM_Is{i} = data.VM_I; 
   % All tday should match
   if i == nFil
       tday = data.tday;
   end
end

% Also get GW data
GWdata = load('VM_GW.mat');
% Should be same k for both datasets
kchoice = GWdata.kchoice;
cd(thisDir);


% Plot both VM_I ratios
figure; cols = {'g', 'b', 'r'};
subplot(2, 1, 1);
hold on;
for i = 1:nMeth
    plot(1:GWdata.nday, GWdata.VM{i}, 'Color', cols{i}, 'LineWidth', 2);
end
hold off; grid off; box off;
xlabel('$s$ (generations)', 'FontSize', 18); 
ylabel(['VM$[I_{s}] \, | \, k$ = ' num2str(kchoice)], 'FontSize', 18);
xlim([2 30]);

subplot(2, 1, 2);
plot(tday, VM_Is{1}, 'Color', cols{1}, 'LineWidth', 2);
hold on; 
plot(tday, VM_Is{2}, 'Color', cols{2}, 'LineWidth', 2);
plot(tday, VM_Is{3}, 'Color', cols{3}, 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$s$ (days)', 'FontSize', 18); 
ylabel(['VM$[I_{s}] \, | \, k$ = ' num2str(kchoice)], 'FontSize', 18);
xlim([20 300]);
cd(loadFol);
saveas(gcf, 'VM_Icomp', 'fig');
cd(thisDir);



