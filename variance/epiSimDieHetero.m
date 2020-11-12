% Simulate dying epidemic via renewal model
function [Iday, Lam, Rsamp, tday, Iwarn] = epiSimDieHetero(nday, ts, Rs, scenNo, distvals, k)

% Assumptions and notes
% - includes heterogeneity via gamma R with dispersion k
% - allows selection of serial interval distributions
% - Rs are distinct reprod nums, ts are switch points 
% - option to remove a startup sequence of zeros
% - warns if sequence of consecutive zero incidences

% All scenarios
scenNam = {'control', 'recovery', 'cascade', 'boom-bust'};

% Variable for true R
Rtrue = zeros(1, nday);
% Change-times in real time
tch = ts;

% Reproduction num profiles
switch(scenNo)
    case 1
        % Rapidly controlled epidemic
        Rtrue(1:tch) = Rs(1);
        Rtrue(tch+1:end) = Rs(2);
    case 2
        % Rapid control that recovers
        Rtrue(1:tch(1)) = Rs(1);
        Rtrue(tch(1)+1:tch(2)) = Rs(2);
        Rtrue(tch(2)+1:end) = Rs(3);
    case 3
        % Two stage control
        Rtrue(1:tch(1)) = Rs(1);
        Rtrue(tch(1)+1:tch(2)) = Rs(2);
        Rtrue(tch(2)+1:end) = Rs(3);
    case 4
        % Exponential rise and fall
        trise = 1:tch; tfall = tch+1:nday;
        % Exponential rise to max at tchange
        Rtrue(trise) =  exp(0.02*(1:tch)); Rmax = Rtrue(tch);
        % Exponential decay from max
        Rtrue(tfall) = Rmax*exp(-0.008*(tfall - tch));
end

% Warning about incidence zeros
Iwarn = 0; % will be set conditionally

% Serial distribution over all tday
serial = serialDistrTypes(nday, distvals);
% Single omega controls distribution
Pomega = serial(1/distvals.omega);

% Daily incidence and infectiousness
Iday = zeros(size(nday)); Lam = Iday; Rsamp = Iday;
% Initialise epidemic
Iday(1) = 10;

% Iteratively generate renewal epidemic
for i = 2:nday
    % Total infectiousness
    Lam(i) = sum(Iday(i-1:-1:1).*Pomega(1:i-1));
    
    % Sample from mean R
    Rsamp(i) = gamrnd(k, Rtrue(i)/k);
    % Renewal incidence
    Iday(i) = poissrnd(Rsamp(i)*Lam(i));
end

% Remove start-up 20 days
idz = 20:nday; tday = idz;
% Adjusted vectors - including tday
Iday = Iday(idz); Rsamp = Rsamp(idz); Lam = Lam(idz);

% Remove small epidemics 
if sum(Iday) < 100 
    Iwarn = 1; 
    disp(['Sum is ' num2str(sum(Iday))]);
end


