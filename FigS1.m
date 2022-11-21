% Multi-deme E optimal design with empirical data
clearvars; clc; close all; tic;

% Assumptions and notes
% - adds time slices across June and July to test resurgence
% - includes max of Rj as another metric
% - accounts for different starting times
% - optimal designs achieved by sampling from qR distributions
% - includes an E and D optimal design for consensus R
% - no interactiion among demes and same serial intervals used

% Directory and where saving (or loading)
thisDir = cd; saveFol = 'Israel/'; 
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
% Extract dates and new cases
Idata = cell(1, nDeme);
for i = 1:nDeme
    Idata{i} = readtable(files(i).name);
end
cd(thisDir);

%% Format empirical data and obtain total infectiousness

% Extract incidence and dates 
nday = length(Idata{1}.date); tday = 1:nday;
Ideme = zeros(nDeme, nday); tdate = Idata{1}.date;
for i = 1:nDeme
    % Assume all files of same length nday
    Ideme(i, :) = Idata{i}.new_cases;
    % Add smoothing (trailing)
    Ideme(i, :) = round(movmean(Ideme(i, :), [6 0]));
end

% Truncate time series so starts from first non-zero term
ndays = zeros(1, nDeme); istarts = zeros(1, nDeme);
for i = 1:nDeme
    istarts(i) = find(Ideme(i, :) > 0, 1, 'first');
    ndays(i) = nday - istarts(i) + 1;
end

% Define a standard serial interval distribution
wmean = 4.7; wvar = 2.9^2;
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

%% Proportion of cases due to Delta variant

% Growth of Delta variant data
cd(saveFol);
vdata = readtable('israel_variant.csv'); 
cd(thisDir);
% Delta and total incidence
Idelta = vdata.Delta; caseDelta = vdata.total;
% Conver to proportion
Idelta = Idelta./caseDelta; tdelta = vdata.Day;

% Extract relevant part within tdate
ndelta = length(Idelta); idelta = zeros(1, ndelta);
for i = 1:ndelta
    idelta(i) = ismember(tdelta(i), tdate);
end
Idelta = Idelta(idelta == 1); tdelta = tdelta(idelta == 1);

% Remove repetitions 
Idelta = Idelta(1:2:end); tdelta = tdelta(1:2:end);
% Get date indices
del1 = find(tdate == tdelta(1)); del2 = find(tdate == tdelta(end));
tdelday = tdate(del1:del2); tdel = tday(del1:del2);

% Duplicate delta cases to get daily rate
ldel = length(tdelday); Idelday = zeros(1, ldel);
for i = 1:ldel
    if ismember(tdelday(i), tdelta)
        % Exact time included
        idd = find(tdelday(i) == tdelta);
        Idelday(i) = Idelta(idd);
    else
        % Keep at last value
        Idelday(i) = Idelday(i-1);
    end
end

% Limit to 95%
tdelday = Idelday(Idelday <= 0.95); tdel = tday(Idelday <= 0.95);
Idelday = Idelday(Idelday <= 0.95);

%% Publishable figure with max of Rj instead of D

% Peak sizes of epidemics
figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
subplot(4, 1, 1);
plot(tdate', Ideme', 'LineWidth', 2);
hold on;
plot(tdate', Itot', 'k', 'LineWidth', 2);
grid off; box off; hold off;
h = gca; h.YScale = 'log';
ylabel('$I_j(t)$', 'FontSize', fnt);
xlim([tdate(8) tdate(end)+1]);

% Get ticks for middle panel
xt = h.XTick; xtlab = h.XTickLabel;
% Find time points of ticks 
tt = zeros(size(xt));
for i = 1:length(xt)-1
    tt(i) = find(tdate == xt(i));
end
tt(length(xt)) = tday(end)+1;

% Key times of interest
tvac1 = xt(4)-2; tvac2 = xt(5)-3;
ttvac1 = tt(4)-2; ttvac2 = tt(5)-3;
h = gca; hold on;
plot([tvac1 tvac1], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
plot([tvac2 tvac2], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
hold off;

subplot(4, 1, [2 3]);
hold on;
plotCIRaw(tday', Rmmaxj', Rlmaxj', Rhmaxj', 'c');
plotCIRaw(tday', RLam', RLaml', RLamh', 'b');
plotCIRaw(tday', RmD', RlD', RhD', 'g');
plotCIRaw(tday', RmE', RlE', RhE', 'r');
plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);
% Addition of variant proportion
stairs(tdel, 3*Idelday, 'LineWidth', 2);
% Vaccination period
h = gca; 
plot([ttvac1 ttvac1], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
plot([ttvac2 ttvac2], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
grid off; box off; hold off;
h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
h.XTick = tt; h.XTickLabel = xtlab;
legend('', '$\max \hat{R}_j(t)$', '', '$\hat{D}(t)$', '', '$\hat{R}(t)$',...
    '', '$\hat{E}(t)$');
legend('boxoff');
ylabel('$\hat{X}(t)$', 'FontSize', fnt);
xlim([tday(8) tday(end)+1]); ylim([0 3]);

subplot(4, 1, 4);
hold on; 
plot(tdate, p1Rjmax, '-', 'color', 'c', 'linewidth', 2);
plot(tdate, p1R, '-', 'color', 'b', 'LineWidth', 2);
plot(tdate, p1D, '-', 'color', 'g', 'linewidth', 2);
plot(tdate, p1E, '-', 'color', 'r', 'LineWidth', 2);
h = gca;
plot([tvac1 tvac1], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
plot([tvac2 tvac2], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
plot(tdate, 0.5*ones(1, nday), '--', 'color', 'k', 'LineWidth', 1);
h = gca; h.YColor = h.XColor;
grid off; box off; hold off;
xlim([tdate(8) tdate(end)+1]); ylim([0 1.05]);
ylabel('P($\hat{X}(t) > 1$)', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt);

if saveFig
    cd(saveFol);
    saveas(gcf, 'Fig5max', 'fig');
    cd(thisDir);
end

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


%% Run 3 slices to test prospective performance

% Define slices across June and July
slices = {'2021-06-10', '2021-06-20', '2021-06-30', '2021-08-01'}; nSlice = length(slices);
% Slice estimate variables


figure('Renderer', 'painters', 'Position', [10 10 800 1000]);
% For each slice truncate time series and estimate
for j = 1:nSlice

    % ID in time series to truncate 
    ndaySlice = find(tdate == datetime(slices{j})); tdaySlice = 1:ndaySlice;
    % Truncate demes and remove those that have not started
    IdemeSlice = Ideme(istarts < ndaySlice, tdaySlice); 
    LdemeSlice = Ldeme(istarts < ndaySlice, tdaySlice);
    nDemeSlice = size(IdemeSlice, 1); iSlice = istarts(istarts < ndaySlice);
    ndaysSlice = ndaySlice-istarts(istarts < ndaySlice)+1;
    ItotSlice = sum(IdemeSlice, 1); LtotSlice = sum(LdemeSlice, 1);

    [Rmmaxj, Rlmaxj, Rhmaxj, RLam, RLaml, RLamh, RmD, RlD, RhD, RmE, RlE, RhE,...
        p1Rjmax, p1R, p1D, p1E]= timeSlice(nDemeSlice, LdemeSlice, IdemeSlice,...
        LtotSlice, ItotSlice, ndaySlice, ndaysSlice, iSlice);

    % Plot of estimates up to slice time
    subplot(nSlice, 2, 2*j-1);
    hold on;
    plotCIRaw(tdaySlice', RLam', RLaml', RLamh', 'b');
    plotCIRaw(tdaySlice', RmE', RlE', RhE', 'r');
    plot(tdaySlice, ones(1, ndaySlice), '--', 'Color', 'k', 'LineWidth', 1);
    % Addition of variant proportion
    stairs(tdel, 3*Idelday, 'k-.', 'LineWidth', 2);
    % Vaccination period
    h = gca;
    plot([ttvac1 ttvac1], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
    grid off; box off; hold off;
    h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
    h.XTick = tt; h.XTickLabel = xtlab;
    ylabel('$\hat{X}(t)$', 'FontSize', fnt);
    xlim([tdaySlice(8) tdaySlice(end)+1]); ylim([0 3]);

    % Plot of p(R > 1) type metrice to slice time
    subplot(nSlice, 2, 2*j);
    hold on;
    plot(tdaySlice', p1R', 'b', 'LineWidth', 2);
    plot(tdaySlice', p1E', 'r', 'LineWidth', 2);
    plot(tdaySlice, 0.5*ones(1, ndaySlice), '--', 'Color', 'k', 'LineWidth', 1);
    % Addition of variant proportion
    stairs(tdel, Idelday, 'k-.', 'LineWidth', 2);
    % Vaccination period
    h = gca;
    plot([ttvac1 ttvac1], h.YLim, '--', 'Color', grey1, 'LineWidth', 1);
    grid off; box off; hold off;
    h = gca; h.YScale = 'linear'; h.YColor = h.XColor;
    h.XTick = tt; h.XTickLabel = xtlab;
    if j == nSlice
        legend('$\hat{R}(t)$', '$\hat{E}(t)$');
        legend('boxoff');
    end
    ylabel('P($\hat{X}(t) > 1$)', 'FontSize', fnt);
    xlim([tdaySlice(8) tdaySlice(end)+1]); ylim([0 1.1]);

end


% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);
