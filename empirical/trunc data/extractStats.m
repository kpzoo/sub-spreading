% Snippet to get means and variances from trunc data
clearvars; clc; close all;

% Assumptions and notes
% - must be in folder with the data

% All data files
filName = dir('sim_*'); nFil = length(filName);
% Means and variances
Imeans = cell(1, nFil); Ivars = Imeans; Rmeans = Imeans; Rvars = Imeans;

% For every data file extract means and vaiances
for i = 1:nFil
    % Main data
    data = load(filName(i).name, 'Iday', 'Rsamp');
    if i == 1
        % Same for all files
        tday = load(filName(i).name, 'tday');
        tday = tday.tday';
    end
    % Statistics
    Imean = mean(data.Iday); Ivar = var(data.Iday);
    Rmean = mean(data.Rsamp); Rvar = var(data.Rsamp);
    % Store in cells to save
    Imeans{i} = Imean; Ivars{i} = Ivar;
    Rmeans{i} = Rmean; Rvars{i} = Rvar;
end

% Save all data
save('statsTraj.mat', 'tday', 'Imeans', 'Ivars', 'Rmeans', 'Rvars', 'filName');