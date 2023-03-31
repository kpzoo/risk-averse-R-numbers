% Multi-deme E optimal design with empirical data
clearvars; clc; close all; tic;

% Assumptions and notes
% - includes max of Rj as another metric
% - accounts for different starting times
% - optimal designs achieved by sampling from qR distributions
% - includes an E and D optimal design for consensus R

% Directory and where saving (or loading)
thisDir = cd; dataNam = 2;
% Decide on datasets to process
switch dataNam
    case 0
        saveFol = 'Case studies/NZ';
        outoforder = 1;
    case 1
        saveFol = 'Case studies/NY';
        outoforder = 1;
    case 2
        saveFol = 'Case studies/Norway';
        outoforder = 0;
    case 3
        saveFol = 'Case studies/Illinois';
        outoforder = 0;
end
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

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

% Remove zero case cities
if dataNam == 2
    nDeme = nDeme - 1;
    % Replace with last one and truncate
    Idata{7} = Idata{end};
end

%% Format empirical data and obtain total infectiousness

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
if outoforder
    tdate = tdate(end:-1:1); Ideme = Ideme(:, end:-1:1);
end

% Truncate time series so starts from first non-zero term
ndays = zeros(1, nDeme); istarts = zeros(1, nDeme);
for i = 1:nDeme
    istarts(i) = find(Ideme(i, :) > 0, 1, 'first');
    ndays(i) = nday - istarts(i) + 1;
end

% Define a standard serial interval distribution
wmean = 4.7; wvar = 2.9^2;
%wmean = 6.5; wvar = wmean*(1/0.65)^2;

% Compose as a gamma distribution
scalePm = wvar/wmean; shapePm = wmean/scalePm;
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

%% EpiFilter estimates of R, D and E numbers

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
clearvars('pstate', 'pstateL', 'pRL', 'pRupL', 'pR', 'pRup');

% Basic D and E optimal design means (no CIs)
Re_D = mean(Rm); Re_E = mean(Rm.^2)./Re_D;
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

%% Publishable figure with max of Rj and D

% Add dates for Norway
if dataNam == 2
    tnpidate = {'2021-09-25', '2021-10-29', '2021-11-15'};
    tnpidate = datetime(tnpidate); 
    ldate = length(tnpidate); tid = zeros(1, ldate);
    for i = 1:ldate
        tid(i) = find(tdate == tnpidate(i));
    end
    % Date for truncation should be 1 Sept
    ttrunc = find(tdate == datetime('2021-09-01'));
else
    ttrunc = 8;
end

% Peak sizes of epidemics
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
subplot(4, 1, 1);
semilogy(tdate', Ideme', 'LineWidth', 2);
hold on;
plot(tdate', Itot', 'k', 'LineWidth', 2);

% Add any npi dates
if dataNam == 2
    for i = 1:length(tnpidate)
        h = gca;
        plot(tdate([tid(i) tid(i)]), h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
    end
end
grid off; box off; hold off;
ylabel('$I_j(t)$', 'FontSize', fnt);
xlim([tdate(ttrunc) tdate(end)+1]);

% Get ticks for middle panel
xt = h.XTick; xtlab = h.XTickLabel;
% Find time points of ticks
tt = zeros(size(xt));
for i = 1:length(xt)-1
    tt(i) = find(tdate == xt(i));
end
tt(length(xt)) = tday(end)+1;

subplot(4, 1, [2 3]);
hold on;
plotCIRaw(tday', Rmmaxj', Rlmaxj', Rhmaxj', 'c');
plotCIRaw(tday', RLam', RLaml', RLamh', 'b');
plotCIRaw(tday', RmD', RlD', RhD', 'g');
plotCIRaw(tday', RmE', RlE', RhE', 'r');
plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);

% Add any npi dates
if dataNam == 2
    for i = 1:length(tnpidate)
        h = gca;
        plot(tday([tid(i) tid(i)]), h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
    end
end

grid off; box off; hold off;
h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
h.XTick = tt; h.XTickLabel = xtlab;
legend('', '$\max \hat{R}_j(t)$', '', '$\hat{R}(t)$', '', '$\hat{D}(t)$',...
    '', '$\hat{E}(t)$');
legend('boxoff');
ylabel('$\hat{X}(t)$', 'FontSize', fnt);
xlim([tday(ttrunc) tday(end)+1]); ylim([0 3]);

subplot(4, 1, 4);
hold on;
plot(tdate, p1Rjmax, '-', 'color', 'c', 'linewidth', 2);
plot(tdate, p1R, '-', 'color', 'b', 'LineWidth', 2);
plot(tdate, p1D, '-', 'color', 'g', 'linewidth', 2);
plot(tdate, p1E, '-', 'color', 'r', 'LineWidth', 2);
plot(tdate, 0.5*ones(1, nday), '--', 'color', 'k', 'LineWidth', 1);

% Add any npi dates
if dataNam == 2
    for i = 1:length(tnpidate)
        plot(tdate([tid(i) tid(i)]), [0 1.05], '--', 'Color', grey1, 'LineWidth', 1);
    end
end

h = gca; h.YColor = h.XColor;
grid off; box off; hold off;
xlim([tdate(ttrunc) tdate(end)+1]); ylim([0 1.05]);
ylabel('P($\hat{X}(t) > 1$)', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt);

if saveFig
    cd(saveFol);
    saveas(gcf, 'Fig5max', 'fig');
    cd(thisDir);
end


% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Remove unneeded variables and save
if saveTrue
    % Variables to save
    varNam = {'tday', 'RmD', 'RlD', 'RhD', 'RLam', 'RLaml', 'RLamh', 'wch', ...
        'RmE', 'RlE', 'RhE', 'Ideme', 'Itot', 'ndays', 'Ltot', 'Ldeme', ...
        'p1R', 'p1D', 'p1E', 'tdate', 'istarts', 'scalePm', 'shapePm'};
    cd(saveFol);
    save(['israel_' num2str(nDeme)  '_data' '.mat'], varNam{:});
    cd(thisDir);
end
