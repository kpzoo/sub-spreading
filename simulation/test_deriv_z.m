% Test variance growth in z with k from NegBin
clearvars; clc; close all; 

% Assumptions and notes
% - tests relationship of var(z) with k and var(z)/E[z]
% - only computes single term derivatives

% Directory of some main code and plotting options
thisDir = cd; cd ..; mainDir = cd; 
mainDir = join([mainDir '/main code']);
cd(thisDir); addpath(mainDir);
% Default plotting options
[grey1, grey2, cmap] = defaultSet(10);

%% Main code for symbolic computation of derivatives

% Smybolic variables
syms('k', 'a', 'real');
% A mean-type function
f = (1 + a/k)^(-k);

% Mean z and its derivative
zm = (1 + a/k)^(-k); dzm = diff(zm);

% Derivative of difference of variance terms
zA = (1 + 2*a/k)^(-k); dzA = diff(zA);
zB = (1 + a/k)^(-2*k); dzB = diff(zB);
dzv = dzA - dzB;

% Dispersion k and constant a
k = logspace(-2, 2, 10000);  
a = 1; lenk = length(k); 

% Examine means 
figure;
semilogx(k, eval(dzm), 'LineWidth', 2);
grid off; box off;
xlabel('$k$'); ylabel('d$\bar{z}$/d$k$');

% Examine variances
figure;
semilogx(k, eval(dzv), 'LineWidth', 2);
grid off; box off;
xlabel('$k$'); ylabel('dV$[z]$/d$k$');

%% Test with a falling a
L = linspace(0, 0.1, 30); L = L(end:-1:1);
L = [L zeros(1, 50)];

% Mean values against k
zmeans = zeros(1, lenk); zvars = zmeans;
zA2 = zmeans; zB2 = zmeans;
for i = 1:lenk
    % Mean terms
    zmean = (1 + L/k(i)).^(-k(i));
    zmeans(i) = prod(zmean);
    zA2 = (1 + 2*L/k(i)).^(-k(i));
    zB2 = (1 + L/k(i)).^(-2*k(i));
    % Var terms
    zvars(i) = prod(zA2) - prod(zB2);
end


figure;
semilogx(k, zmeans, 'LineWidth', 2);
grid off; box off;
xlabel('$k$'); ylabel('$\bar{z}$');


figure;
semilogx(k, zvars, 'LineWidth', 2);
grid off; box off;
xlabel('$k$'); ylabel('V$[z]$');
