% Process case study data into reproduction number estimates
function [RL, RLl, RLh, p1R, RmE, RlE, RhE, p1E, tdate, tday, Ideme, Ldeme,...
    Itot, Ltot] = procEmpCaseStudy(saveFol, caseID, nsamps, scalePm, shapePm, thisDir)

% Assumptions and notes
% - performs both filtering and smoothing and samples
% - outputs the mean and CI of E and R

%% Extracting and cleaning data to perform estimation

% Load incidence data and serial intervals
cd(saveFol);
% Read and count files with incidence data
files = dir('city*'); nDeme = length(files);

% Extract dates and new cases and sum
Idata = cell(1, nDeme); Isums = zeros(1, nDeme);
for i = 1:nDeme
    Idata{i} = readtable(files(i).name);
    Isums(i) = sum(Idata{i}.new_cases);
end
cd(thisDir);

% Remove zero case cities from Norway data
if caseID == 2
    nDeme = nDeme - 1;
    % Replace with last one and truncate
    Idata{7} = Idata{end};
end

% Truncate time series
nday = length(Idata{1}.date); tday = 1:nday;
% Extract incidence and dates
Ideme = zeros(nDeme, nday); tdate = Idata{1}.date(tday);
for i = 1:nDeme
    % Assume all files of same length nday
    Ideme(i, :) = Idata{i}.new_cases(tday);
    % Add smoothing (trailing)
    Ideme(i, :) = round(movmean(Ideme(i, :), [6 0]));
    % Remove any negative values
    Ideme(i, Ideme(i, :) < 0)= 0;
end

% Dates out of order
if ismember(caseID, [3 4])
    tdate = tdate(end:-1:1); Ideme = Ideme(:, end:-1:1);
end

% Truncate time series so starts from first non-zero term
ndays = zeros(1, nDeme); istarts = zeros(1, nDeme);
for i = 1:nDeme
    istarts(i) = find(Ideme(i, :) > 0, 1, 'first');
    ndays(i) = nday - istarts(i) + 1;
end

% Serial interval distribution
wch = gamcdf(tday, shapePm, scalePm) - gamcdf(tday-1, shapePm, scalePm);

% Total infectiousness of each deme
Ldeme = Ideme;
for j = 1:nDeme
    for i = 2:nday
        Ldeme(j, i) = sum(Ideme(j, i-1:-1:1).*wch(1:i-1));
    end
end

% Aggregate the epidemics and also Lam as serial interval fixed
Itot = sum(Ideme, 1); Ltot = sum(Ldeme, 1);

%% EpiFilter estimates of R and E numbers

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1; m = 1000; 
% Uniform prior over grid of size m
p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% Estimates variables for EpiFilter
Rm = zeros(nDeme, nday); Rl = Rm; Rh = Rm; 
% Smoothed and filtered distributions by deme
qR = cell(1, nDeme); pR = qR; 

% Smoothed estimates and distributions
for i = 1:nDeme
    % EpiFilter estimate from each deme
    [~, ~, ~, ~, pR{i}, pRup, pstate] = runEpiFilter(Rgrid, m, eta, ndays(i), p0,...
        Ldeme(i, istarts(i):end), Ideme(i, istarts(i):end));
    [~, Rl(i, istarts(i):end), Rh(i, istarts(i):end), Rm(i, istarts(i):end),...
        qR{i}] = runEpiSmoother(Rgrid, m, ndays(i), pR{i}, pRup, pstate);

    disp(['Completed deme ' num2str(i) ' of ' num2str(nDeme)]);
end

% E optimal distribution by sampling 
RmE = zeros(2, nday); RlE = RmE; RhE = RmE; 
RL = RmE; RLl = RmE; RLh = RmE; 
p1E = RmE; p1R = RmE;

% Filtered and smoothed overall E estimates
[RmE(1,:), RlE(1,:), RhE(1,:), p1E(1,:), ~] = sampEDist(nsamps, nday, pR, nDeme, istarts, Rgrid);
[RmE(2,:), RlE(2,:), RhE(2,:), p1E(2,:), ~] = sampEDist(nsamps, nday, qR, nDeme, istarts, Rgrid);

% Aggregrate R estimates over demes
[~, RLl(1,:), RLh(1,:), RL(1,:), pRL, pRupL, pstateL] = runEpiFilter(Rgrid,...
    m, eta, nday, p0, Ltot, Itot);
[~, RLl(2,:), RLh(2,:), RL(2,:), qRL] = runEpiSmoother(Rgrid, m, nday, pRL, pRupL, pstateL);

% Prob R > 1
id1 = find(Rgrid <= 1, 1, 'last');
for i = 1:nday
    p1R(1, i) = 1 - sum(pRL(i, 1:id1), 2)';
    p1R(2, i) = 1 - sum(qRL(i, 1:id1), 2)';
end