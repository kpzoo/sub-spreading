% Get truncation point in gamma to satisfy a mean 
function [trDist, trStat] = optGamTruncFix(truncType, a, b, rhoMean, Rmean, thresh)

% Assumptions and notes
% - threshold fixed so truncated a, b modified
% - [a, b] must be k and Rmean*rhomean/k
% - undersampling or control causes a truncation of the R-gamma distrib
% - find the theshold value such that a desired target mean R is obtained

switch(truncType)
    case 1
        % Lower end truncated due to undersampling
         
        % Start with original distribution
        pd = makedist('Gamma', 'a', a, 'b', b);
        % For fixed threshold find scaling to obtain target mean
        targMean = rhoMean*Rmean;
        
        % Aim to truncate so mean goes up to target
        low = @(x) lowerTrunc(x, pd, targMean, thresh);
        % Optimum parameter
        options = optimset('TolFun',1e-5);
        [scale, fval] = fminsearch(low, 0.1, options);
        % Confirm result 
        if fval > 10^-3
            assignin('base', 'fval', fval);
            assignin('base', 'scale', scale);
            error('Issue with minimisation');
        else
            % Output truncated distribution
            pd.b = pd.b/scale;
            trDist = truncate(pd, thresh, inf);
            trStat.mean = mean(trDist);
            trStat.var = var(trDist);
        end
        
    case 2
        % Upper end truncated due to control

end



%% Lower truncation function
function err = lowerTrunc(x, pd, targMean, thresh)

% Re-scale distribution
pd.b = pd.b/x;

% Truncated distribution and mean
tr = truncate(pd, thresh, inf);
gmean = mean(tr);

% Squared error of means
err = 100*abs(1 - gmean/targMean);

