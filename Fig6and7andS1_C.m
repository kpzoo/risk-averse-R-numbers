% Empirical case studies with filtering and crossing times
clearvars; clc; close all; tic;

% Assumptions and notes
% - 6 case studies and only smoothed E and R shown
% - then zoom into time slices for resurgence
% - accounts for different starting times
% - optimal designs achieved by sampling from qR distributions

% Directory and where saving (or loading)
thisDir = cd; saveFol = 'Case studies/';
% Booleans for saving 
saveTrue = 0; saveFig = 0;

% All case study names
dataNam = {'Israel', 'Norway', 'NZ', 'NY', 'Illinois', 'UK region'};
labNam = {'Israel', 'Norway', 'New Zealand', 'New York', 'Illinois', 'UK regions'};
% Num of studies and samples
nStudy = length(dataNam); nsamps = 50000;
 
% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% Define a standard serial interval distribution
wmean = 4.7; wvar = 2.9^2;
% Compose as a gamma distribution
scalePm = wvar/wmean; shapePm = wmean/scalePm;


%% Obtain filtered and smoothed estimates for each dataset

% Store results from datasets
RL = cell(1, nStudy); RLl = RL; RLh = RL; p1R = RL; RmE = RL; RlE = RL; RhE = RL;
p1E = RL; tdate = RL; tday = RL; Ideme = RL; Ldeme = RL; Itot = RL; Ltot = RL;

% Posterior estimates with 95% CIs
for ii = 1:nStudy
    saveFolii = [saveFol dataNam{ii}];

    % Main function to extract data and estimates
    [RL{ii}, RLl{ii}, RLh{ii}, p1R{ii}, RmE{ii}, RlE{ii}, RhE{ii}, p1E{ii},...
        tdate{ii}, tday{ii}, Ideme{ii}, Ldeme{ii}, Itot{ii}, Ltot{ii}]...
        = procEmpCaseStudy(saveFolii, ii, nsamps, scalePm, shapePm, thisDir);

    disp(['Completed ' num2str(ii) ' of ' num2str(nStudy)]);
end

% Times when lower interval estimate crosses 1 and Prob > 1 crosses 0.95
tlowR = cell(1, nStudy); tlowE = tlowR; t1R = tlowR; t1E = tlowR;
for ii = 1:nStudy
    [tlowR{ii}, t1R{ii}] = getChgPt(tday{ii}, RLl{ii}(1, :), p1R{ii}(1, :));
    [tlowE{ii}, t1E{ii}] = getChgPt(tday{ii}, RlE{ii}(1, :), p1E{ii}(1, :));
end

% Smoothed times when lower interval estimate crosses 1 and Prob > 1 crosses 0.95
tlowRs = cell(1, nStudy); tlowEs = tlowRs; t1Rs = tlowRs; t1Es = tlowRs;
for ii = 1:nStudy
    [tlowRs{ii}, t1Rs{ii}] = getChgPt(tday{ii}, RLl{ii}(2, :), p1R{ii}(2, :));
    [tlowEs{ii}, t1Es{ii}] = getChgPt(tday{ii}, RlE{ii}(2, :), p1E{ii}(2, :));
end

% Define times for each case study where to investigate resurgence
tstartdate = {'2021-06-08', '2021-07-15', '2022-01-04', '2020-09-03', '2021-06-21', '2021-05-16'};
tstartdate = datetime(tstartdate); tenddate = tstartdate + 15;
% Find ids in date time series for endpoints
idstart = zeros(1, nStudy); idend = idstart;
for ii = 1:nStudy
    idstart(ii) = find(tdate{ii} == tstartdate(ii));
    idend(ii) = find(tdate{ii} == tenddate(ii));
end

% Times for Israel
t1Is = zeros(2, 2); R1Is = t1Is; thr1 = [0.5 0.95]; tid1 = idstart(1):idend(1); 
for i = 1:2
    t1Is(1, i) = find(p1R{1}(1,tid1) >= thr1(i), 1, "first") + tid1(1) - 1;
    t1Is(2, i) = find(p1E{1}(1,tid1) >= thr1(i), 1, "first") + tid1(1) - 1;
end
R1Is(1, 1) = find(RLl{1}(1,tid1) >= 1, 1, "first") + tid1(1) - 1;
R1Is(1, 2) = find(RL{1}(1,tid1) >= 1, 1, "first") + tid1(1) - 1;
R1Is(2, 1) = find(RlE{1}(1,tid1) >= 1, 1, "first") + tid1(1) - 1;
R1Is(2, 2) = find(RmE{1}(1,tid1) >= 1, 1, "first") + tid1(1) - 1;
% Convert to dates
t1Is = tdate{1}(t1Is); R1Is = tdate{1}(R1Is);

%% Publishable figure with all estimates and change times

% Simple figure of cases by demes for supplement
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
for ii = 1:nStudy
    subplot(round(nStudy)/2, 2, ii);
    yyaxis('left'); hold on;
    % Incidence for all demes and n
    plot(tdate{ii}, Ideme{ii}, '-', 'Color', grey1, 'LineWidth', 2);
    grid off; box off; hold off;
    xlim([tdate{ii}(8) tdate{ii}(end)+1]);
    h = gca; h.YAxis(1).TickLabelColor = h.XAxis.TickLabelColor;
    ylabel(['$I_j(t) | p = $ ' num2str(size(Ideme{ii}, 1))], 'FontSize', fnt);
    h.YLabel.Color = h.XAxis.TickLabelColor;
    title(labNam{ii}, FontSize=fnt);

    yyaxis('right'); hold on;
    plot(tdate{ii}, Itot{ii}, 'k-', 'LineWidth', 2);
    plot([tstartdate(ii) tstartdate(ii)], h.YLim, '--', 'Color', grey1, 'LineWidth', 2);
    plot([tenddate(ii) tenddate(ii)], h.YLim, '--', 'Color', grey1, 'LineWidth', 2);
    h = gca; h.YAxis(2).TickLabelColor = h.XAxis.TickLabelColor;
    h.YAxis(2).TickLabels = ''; grid off; box off; hold off;
end
if saveFig
    cd(saveFol);
    saveas(gcf, 'FigS1_A', 'fig');
    saveas(gcf, 'FigS1_A', 'png');
    cd(thisDir);
end

% Full smoothed estimates of all case studies
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
for ii = 1:nStudy
    subplot(round(nStudy)/2, 2, ii);
    % R and E estimates
    yyaxis('left'); hold on;
    plot(tdate{ii}, ones(1, length(tday{ii})), '--', 'Color', 'k', 'LineWidth', 1);
    plotCIRaw(tday{ii}', RL{ii}(2,:)', RLl{ii}(2,:)', RLh{ii}(2,:)', 'b');
    plotCIRaw(tday{ii}', RmE{ii}(2,:)', RlE{ii}(2,:)', RhE{ii}(2,:)', 'r');
    grid off; box off; hold off;
    xlim([tdate{ii}(8) tdate{ii}(end)+1]);
    h = gca; h.YAxis(1).TickLabelColor = h.XAxis.TickLabelColor;
    % Label and title by country
    h.YLabel.Color = h.XAxis.TickLabelColor;
    ylabel('$\hat{X}(t)$', 'FontSize', fnt);
    title(labNam{ii}, FontSize=fnt);

    % Total incidence
    yyaxis('right'); hold on; h = gca;
    plot(tdate{ii}, Itot{ii}, 'k-', 'LineWidth', 2);
    plot([tstartdate(ii) tstartdate(ii)], h.YLim, '--', 'Color', grey1, 'LineWidth', 2);
    plot([tenddate(ii) tenddate(ii)], h.YLim, '--', 'Color', grey1, 'LineWidth', 2);
    h = gca; h.YAxis(2).TickLabelColor = h.XAxis.TickLabelColor;
    h.YAxis(2).TickLabels = ''; grid off; box off; hold off;
end
if saveFig
    cd(saveFol);
    saveas(gcf, 'Fig6', 'fig');
    saveas(gcf, 'Fig6', 'png');
    cd(thisDir);
end

% Plot probabilities > 1 over regions of interest
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
for ii = 1:nStudy
    subplot(round(nStudy)/2, 2, ii);
    tid = idstart(ii):idend(ii);

    % Get intersects with thresholds
    tintR = find(p1R{ii}(1,tid)-0.95 == min(abs(p1R{ii}(1,tid) - 0.95)));

    % P(R > 1) and P(E > 1)
    yyaxis('left'); hold on;
    plot(tdate{ii}(tid), p1R{ii}(1,tid), '--', 'color', 'b', 'LineWidth', 2);
    plot(tdate{ii}(tid), p1E{ii}(1,tid), '--', 'color', 'r', 'LineWidth', 2);
    plot(tdate{ii}(tid), 0.95*ones(1, length(tid)), '--', 'color', 'k', 'LineWidth', 1);
    h = gca; h.YAxis(1).TickLabelColor = h.XAxis.TickLabelColor;
    % Label and title by country
    h.YLabel.Color = h.XAxis.TickLabelColor;
    ylabel('P$(\hat{X}(t) > 1)$', 'FontSize', fnt);
    title(labNam{ii}, FontSize=fnt);
    ylim([h.YLim(1), 1.1]);

    % R and E estimates at lower CI
    yyaxis('right'); hold on; 
    plot(tdate{ii}(tid), RLl{ii}(1,tid), '-', 'color', 'b', 'LineWidth', 2);
    plot(tdate{ii}(tid), RlE{ii}(1,tid), '-', 'color', 'r', 'LineWidth', 2);
    plot(tdate{ii}(tid), ones(1, length(tid)), '-', 'color', 'k', 'LineWidth', 1);
    h = gca; h.YAxis(2).TickLabelColor = h.XAxis.TickLabelColor;
    ylabel('$\hat{X}(t)$', 'FontSize', fnt);
    h.YLabel.Color = h.XAxis.TickLabelColor;
    grid off; box off; hold off;
end
if saveFig
    cd(saveFol);
    saveas(gcf, 'Fig7', 'fig');
    saveas(gcf, 'Fig7', 'png');
    cd(thisDir);
end

% Timing and data saving
tsim = toc/60;
disp(['Run time = ' num2str(tsim)]);

% Remove unneeded variables and save
if saveTrue
    % Variables to save
    varNam = {'tday', 'RL', 'RLl', 'RLh', 'RmE', 'RlE', 'RhE', 'Ideme', 'Itot',...
        'ndays', 'Ltot', 'Ldeme','p1R', 'p1E', 'tdate', 'scalePm', 'shapePm', 'wch'};
    cd(saveFol);
    save(['comp_' num2str(nStudy) '.mat'], varNam{:});
    cd(thisDir);
end
