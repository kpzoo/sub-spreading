% Snippet to get VM ratios from trunc data after extractStats
clearvars; clc; close all;

% Assumptions and notes
% - must be in folder with the data statsTraj

% Possibe controls
ctrlName = {'uniform', 'sub-spreading', 'super-spreading'};
ctrlID = {'unif', 'sub', 'sup'}; nCtrl = length(ctrlID);
cols = {'g', 'b', 'r'};

% All data files
fileName = dir('statsTraj_*'); nFil = length(fileName);

% Process in order of control names
ord = zeros(1, nFil);
for i = 1:nFil
    nam = fileName(i).name;
    % Find which id it belongs to
    id = zeros(1, nCtrl);
    for j = 1:nCtrl
        id(j) = ~isempty(regexp(nam, ctrlID{j}, 'ONCE'));
    end
    % Assign order
    ord(i) = find(id);
end

% Reorder files
fileTemp = fileName;
for i = 1:nFil
    fileName(ord(i)) = fileTemp(i);
end

% Main VM ratios to consider
VM_I = cell(1, 1); VM_R = VM_I;

% Load data from each file 
for i = 1:nFil
    data = load(fileName(i).name, 'Imeans', 'Ivars', 'Rmeans', 'Rvars');
    % Each variable is for several k
    VM_Itemp = cell(1, length(data.Imeans)); VM_Rtemp = VM_Itemp;
    for j = 1:length(data.Imeans)
        VM_Itemp{j} = data.Ivars{j}./data.Imeans{j};
        VM_Rtemp{j} = data.Rvars{j}./data.Rmeans{j};
    end
    % Store all results
    VM_I{i} = VM_Itemp; VM_R{i} = VM_Rtemp;
    clear VM_Itemp VM_Rtemp
end

% Plot VM ratios for a given k, kid
kid = 12;

figure;
subplot(2, 1, 1);
semilogy(VM_I{1}{kid}, 'Color', cols{1}, 'LineWidth', 2);
hold on;
semilogy(VM_I{3}{kid}, 'Color', cols{3}, 'LineWidth', 2);
semilogy(VM_I{2}{kid}, 'Color', cols{2}, 'LineWidth', 2);
grid off; box off; hold off;
xlabel('$s$ (days)', 'FontSize', 18);
ylabel('VM[$I_s$]', 'FontSize', 18);
subplot(2, 1, 2);
semilogy(VM_R{1}{kid}, 'Color', cols{1}, 'LineWidth', 2);
hold on;
semilogy(VM_R{3}{kid}, 'Color', cols{3}, 'LineWidth', 2);
semilogy(VM_R{2}{kid}, 'Color', cols{2}, 'LineWidth', 2);
grid off; box off; hold off;
xlabel('$s$ (days)', 'FontSize', 18);
ylabel('VM[$R_s$]', 'FontSize', 18);
