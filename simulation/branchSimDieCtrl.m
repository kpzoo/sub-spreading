% Simulate dying epidemic via Galton Watson branching process
function [Iday, Iwarn] = branchSimDieCtrl(nday, ts, trDist1, trDist2)

% Assumptions and notes
% - includes heterogeneity via gamma R with dispersion k
% - applies different styles of control to this dispersion
% - only does a step change in reproduction numbers

% Daily incidence and infectiousness
Iday = zeros(1, nday); 
% Initialise epidemic and warning of 0s
Iday(1) = 10; Iwarn = 0;

% Iteratively generate branching process
i = 2;
while(Iday(i-1) > 0 && i <= nday)
    % Draw offspring samples
    if i <= ts(1)
        Iday(i) = poissrnd(sum(random(trDist1, [1 Iday(i-1)])));
    else
        Iday(i) = poissrnd(sum(random(trDist2, [1 Iday(i-1)])));
    end
    i = i + 1;
end

% Remove small epidemics
if sum(Iday) < 100
    Iwarn = 1;
    disp(['Sum is ' num2str(sum(Iday))]);
end


