% True probability that epidemic is eliminated given Rtrue
function [z0, zseq0] = getProbEliminTrueVecFast(I, Rtrue, PomegaMax, zlen)

% Assumptions and notes
% - hopefully faster than getProbEliminTrueVector via PomegaMax
% - R is known, a vector of samples is provided
% - R vector must match Lcurr(ir+1) in size
% - append with pseudo-data of 0s to compute
% - only input data up to point want to check elimination
% - posterior for R based on k days look-back

% Define a zero look-ahead sequence
Iz = zeros(1, zlen); idz = length(I) + 1;

% Append epidemic curve with pseudo-data
Icurr = [I Iz]; ncurr = length(Icurr);

% Gamma serial distribution
Pomega = PomegaMax(1:ncurr);

% % Compute successive total infectiousness for this I
% Lcurr = zeros(1, ncurr);
% for i = 2:ncurr
%     % Relevant part of SI: Pomega(1:i-1))
%     Lcurr(i) = Icurr(i-1:-1:1)*Pomega(1:i-1)';
% end

% Call mex file
Lcurr = succLam_mex(Icurr, Pomega, ncurr);

% Range of index time points
ir = idz:ncurr-1; 

% Probability of elimination sequence
zseq0 = exp(-Rtrue(1:end-1).*Lcurr(ir+1));
% Probability of elimination
z0 = prod(zseq0);





