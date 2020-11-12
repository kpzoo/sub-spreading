% Simulate dying epidemic via Galton Watson branching process
function [Iday, Iwarn] = branchSimDieCtrlTwo(nday, ts, trDist1, trDist2)

% Assumptions and notes
% - includes an extra generation dependence
% - includes heterogeneity via gamma R with dispersion k
% - applies different styles of control to this dispersion
% - only does a step change in reproduction numbers

% Daily incidence and infectiousness
Iday = zeros(1, nday); 
% Initialise epidemic and warning of 0s
Iday(1) = 10; Iwarn = 0;

% Iteratively generate branching process
i = 2;
while(i <= nday)
    % Draw offspring samples
    if i <= ts(1)
        Iday(i) = sum(random(trDist1, [1 Iday(i-1)]));
    else
        Iday(i) = sum(random(trDist2, [1 Iday(i-1)]));
    end
    % Extra dependence
    Iday(i) = poissrnd(0.5*(Iday(i) + Iday(i-1)));
    i = i + 1;
end

% Remove small epidemics
if sum(Iday) < 100
    Iwarn = 1;
    disp(['Sum is ' num2str(sum(Iday))]);
end


