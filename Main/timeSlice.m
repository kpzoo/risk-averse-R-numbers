% Function to perform time slice estimates (R and E)
function [Rmmaxj, Rlmaxj, Rhmaxj, RLam, RLaml, RLamh, RmD, RlD, RhD, RmE, RlE, RhE,...
    p1Rjmax, p1R, p1D, p1E]= timeSlice(nDeme, Ldeme, Ideme, Ltot, Itot, nday, ndays, istarts)

% Assumptions and notes
% - setup for Israel case study but generalisable
% - run EpiFilter then sample distributions up to time slice


% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% Estimates variables for EpiFilter
Rm = zeros(nDeme, nday); Rl = Rm; Rh = Rm; qR = cell(1, nDeme);

% Smoothed estimates and distributions
for i = 1:nDeme
    % EpiFilter estimate from each deme
    [~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, ndays(i), p0,...
        Ldeme(i, istarts(i):end), Ideme(i, istarts(i):end));
    [~, Rl(i, istarts(i):end), Rh(i, istarts(i):end), Rm(i, istarts(i):end),...
        qR{i}] = runEpiSmoother(Rgrid, m, ndays(i), pR, pRup, pstate);
    
    disp(['Completed deme ' num2str(i) ' of ' num2str(nDeme)]);
end

% Aggregrate estimate over demes
[~, ~, ~, ~, pRL, pRupL, pstateL] = runEpiFilter(Rgrid, m, eta, nday, p0, Ltot, Itot);
[~, RLaml, RLamh, RLam, qRLam] = runEpiSmoother(Rgrid, m, nday, pRL, pRupL, pstateL);

% Prob > 1 for metrics
p1E = zeros(1, nday); p1D = p1E; p1R = p1D; p1Rjmax = p1D;

% D and E optimal distribution by sampling and max Rj
nsamps = 50000; RmD = zeros(1, nday); RlD = RmD; RhD = RmD; RmE = RmD; 
RlE = RmD; RhE = RmD; Rhmaxj = RmD; Rlmaxj = RmD; Rmmaxj = RmD;

for i = 1:nday
    % Individual distribution samples
    xDeme = zeros(nDeme, nsamps);
    for j = 1:nDeme
        if i >= istarts(j)
            xDeme(j, :) = datasample(Rgrid, nsamps, 'Weights', qR{j}(i-istarts(j)+1, :));
        end
    end
    % D and E optimal samples for this day
    Dsamp = mean(xDeme); 
    Esamp = mean(xDeme.^2)./Dsamp;
    % Sample maximum across demes
    Rjmaxsamp = max(xDeme);

    % Statistics of D designs
    RmD(i) = mean(Dsamp); 
    Dquants = quantile(Dsamp, [0.025, 0.975]);
    RlD(i) = Dquants(1); RhD(i) = Dquants(2);

     % Statistics of E designs
    RmE(i) = mean(Esamp);
    Equants = quantile(Esamp, [0.025, 0.975]);
    RlE(i) = Equants(1); RhE(i) = Equants(2);

    % Statistics of max of Rj metric
    Rmmaxj(i) = mean(Rjmaxsamp);
    Rjmaxquants = quantile(Rjmaxsamp, [0.025, 0.975]);
    Rlmaxj(i) = Rjmaxquants(1); Rhmaxj(i) = Rjmaxquants(2);

    % Prob of D > 1, E > 1, max Rj > 1
    p1D(i) = length(Dsamp(Dsamp > 1))/nsamps;
    p1E(i) = length(Esamp(Esamp > 1))/nsamps;
    p1Rjmax(i) = length(Rjmaxsamp(Rjmaxsamp > 1))/nsamps;
end

id1 = find(Rgrid <= 1, 1, 'last');
for i = 1:nday
    p1R(i) = 1 - sum(qRLam(i, 1:id1), 2)';
end