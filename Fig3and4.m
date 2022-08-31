% Two deme renewal models and estimates
clearvars; clc; close all; tic;

% Assumptions and notes
% - optimal designs achieved by sampling from qR distributions
% - includes an E and D optimal design for consensus R
% - compare local and aggregrated R estimates and I predictions
% - same serial interval and total time period for epidemics
% - no interaction between demes i.e. travel bans

% Directory and where saving
thisDir = cd; saveFol = 'Results/'; 
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

%% Simulate epidemics for both demes

% Define number of days to simulate
nday0 = 301; tday0 = 1:nday0;
% Defined cases for figures
figNo = 4; distNo = 2;
switch(figNo)
    case 3
        % Step change resurgence
        scenNos = [2 1 2]; 
        offsets = [1 30 60];
    case 4
        % Sinusoidal resurgence example
        scenNos = [6 6 4]; 
        offsets = [1 60 1];
end

% Num demes considered
nDeme = length(scenNos); disp(['Examining ' num2str(nDeme) ' demes']);
% Output incidence curve variables
Ideme = cell(1, nDeme); Ldeme = Ideme; Rtrue = Ideme; 

for i = 1:nDeme
    % Simulation parameters and warning
    Iwarn = 1; simVals = setupScenario(scenNos(i));
    simVals.offset = offsets(i);
    % Simulate epidemic scenarios and truncate
    while Iwarn
        [Ideme{i}, Ldeme{i}, Rtrue{i}, tday, Iwarn, distvals] = ...
            epiSimScenDeme(scenNos(i), nday0, distNo, simVals);
    end
    if Iwarn
        warning('Sequences of zero incidence');
    end
    % Truncated observation period - same for all demes
    nday = length(tday);
end

% Saving data and figs
namstr = [num2str(nday) '_' num2str(nDeme) '_' num2str(distNo)];

% Convert individual epidemics to matrices for ease
Ideme = cell2mat(Ideme'); Ldeme = cell2mat(Ldeme'); 
% True deme R and average
Rtrue = cell2mat(Rtrue'); Ravg = mean(Rtrue, 1);

% Aggregate the epidemics and also Lam as serial interval fixed
Itot = sum(Ideme, 1); Ltot = sum(Ldeme, 1);

%% Estimate transmission among epidemics in demes

% Grid limits and noise level
Rmin = 0.01; Rmax = 10; eta = 0.1;
disp(['[Eta, Rmax, Rmin] = ' num2str(eta) ', ' num2str(Rmin) ', ' num2str(Rmax)]);
% Uniform prior over grid of size m
m = 1000; p0 = (1/m)*ones(1, m); Rgrid = linspace(Rmin, Rmax, m);

% Estimates variables for EpiFilter
Rm = zeros(nDeme, nday); Rl = Rm; Rh = Rm; qR = cell(1, nDeme);

% Smoothed estimates and distributions
for i = 1:nDeme
    % EpiFilter estimate from each deme
    [~, ~, ~, ~, pR, pRup, pstate] = runEpiFilter(Rgrid, m, eta, nday, p0, Ldeme(i, :), Ideme(i, :));
    [~, Rl(i, :), Rh(i, :), Rm(i, :), qR{i}] = runEpiSmoother(Rgrid, m, nday, pR, pRup, pstate);
    
    disp(['Completed deme ' num2str(i) ' of ' num2str(nDeme)]);
end

% Aggregrate estimate over demes
[~, ~, ~, ~, pRL, pRupL, pstateL] = runEpiFilter(Rgrid, m, eta, nday, p0, Ltot, Itot);
[~, RLaml, RLamh, RLam, qRLam] = runEpiSmoother(Rgrid, m, nday, pRL, pRupL, pstateL);

% Basic D and E optimal design means (no CIs)
Re_D = mean(Rm); Re_E = mean(Rm.^2)./Re_D;

% D and E optimal distribution by sampling
nsamps = 10000; RmD = zeros(1, nday); RlD = RmD; RhD = RmD;
RmE = RmD; RlE = RmD; RhE = RmD;
for i = 1:nday
    % Individual distribution samples
    xDeme = zeros(nDeme, nsamps);
    for j = 1:nDeme
        xDeme(j, :) = datasample(Rgrid, nsamps, 'Weights', qR{j}(i, :));
    end
    % D and E optimal samples for this day
    Dsamp = mean(xDeme); Esamp = mean(xDeme.^2)./Dsamp;

    % Statistics of D designs
    RmD(i) = mean(Dsamp); 
    Dquants = quantile(Dsamp, [0.025, 0.975]);
    RlD(i) = Dquants(1); RhD(i) = Dquants(2);

     % Statistics of E designs
    RmE(i) = mean(Esamp);
    Equants = quantile(Esamp, [0.025, 0.975]);
    RlE(i) = Equants(1); RhE(i) = Equants(2);
end

%% Publishable figures

% Axes handles
ax1 = zeros(1, nDeme); ax2 = ax1; grey3 = [0.9 0.9 0.9];
% Comparison of each deme estimates
figure('Renderer', 'painters', 'Position', [10 10 700 1100]);
for i = 1:nDeme
    ax1(i) = subplot(nDeme+2, 1, i);
    yyaxis('right');
    stairs(tday, Ideme(i, :), 'Color',  cmap(4, :), 'LineWidth', 2);
    h = gca; h.YColor = h.XColor;
    yyaxis('left');
    plotCIRaw(tday', Rm(i, :)', Rl(i, :)', Rh(i, :)', 'b');
    hold on; h = gca; h.YColor = h.XColor;
    stairs(tday, Rtrue(i, :), '--', 'Color', 'k', 'linewidth', 2);
    plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);
    grid off; box off; hold off; xlim([tday(2) tday(end)]);
    ylabel(['$\hat{R}_{' num2str(i) '}(t)$'], 'FontSize', fnt);
end
% Total and consensus R estimates
ax1(nDeme+1) = subplot(nDeme+2, 1, [nDeme+1 nDeme+2]);
yyaxis('right');
stairs(tday, Itot, 'Color', cmap(4, :), 'LineWidth', 2);
h = gca; h.YColor = h.XColor;
yyaxis('left');
plotCIRaw(tday', RmD', RlD', RhD', 'g');
hold on; h = gca; h.YColor = h.XColor;
plotCIRaw(tday', RLam', RLaml', RLamh', 'b');
plotCIRaw(tday', RmE', RlE', RhE', 'r');
plot(tday, ones(1, nday), '--', 'Color', 'k', 'LineWidth', 1);
grid off; box off; hold off;
ylabel('$\hat{R}(t), \hat{D}(t), \hat{E}(t)$', 'FontSize', fnt);
xlabel('$t$ (days)', 'FontSize', fnt);
xlim([tday(2) tday(end)]); linkaxes(ax1, 'y'); 

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Remove unneeded variables and save
if saveTrue
    cd(saveFol);
    save(['deme' namstr  'data' '.mat']);
    cd(thisDir);
end
