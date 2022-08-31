% Simple tests of the sensitivity of D, R and E
clearvars; clc; close all; tic;

% Assumptions and notes
% - optimal designs achieved by sampling EpiEstim gamma posteriors
% - examine at two scenarios of L1 = 20 and 100

% Directory and where saving
thisDir = cd; saveFol = 'Results/';
% Booleans for saving
saveTrue = 0; saveFig = 0;

% Directory of some main code and plotting options
cd('Main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);


%% Show E more sensitive than R but less sensitive than R1

% Examine several group 2 active infection sizes
aj = [0.5 1 2]; lena = length(aj); ndel = 20;
% Delta to incidence on region 1 and sample size
del = linspace(0, 1, ndel); nSamp = 20000;


% Estimates at different L1 values
L1 = [20 100]; len1 = length(L1); ms = cell(1, len1); vs = ms; Fs = ms;
for j = 1:len1
    % Mean, var, P>1 and domain x
    [ms{j}, vs{j}, Fs{j}, x] = getERgamma(aj, L1(j), ndel, lena);
end

% Construct plot from estimates directly
figure('Renderer', 'painters', 'Position', [10 10 800 800]);
for i = 1:lena
    % Plots of relative prob resurgence
    subplot(3, 2, 1+2*(i-1));
    % Prob resurgence
    yyaxis('left');
    hold on; j = 1;
    %plot(x(i, :), Fs{j}.R(i, :)./Fs{j}.R1(i, :) - 1, 'b-', 'LineWidth', 2);
    %plot(x(i, :), Fs{j}.E(i, :)./Fs{j}.R1(i, :) - 1, 'r-', 'LineWidth', 2);
    plot(x(i, :), Fs{j}.R1(i, :), 'k-', 'LineWidth', 2);
    plot(x(i, :), Fs{j}.R(i, :), 'b-', 'LineWidth', 2);
    plot(x(i, :), Fs{j}.E(i, :), 'r-', 'LineWidth', 2);
    hold off; box off; grid off;
    h = gca; h.YColor = h.XColor; 
    ymin = h.YLim(1); ylim([ymin 1.01]);
    ylabel(['P$(X > 1) | r = $ ' num2str(aj(i))], 'FontSize', fnt);
    if i == 1
        title(['$\Lambda_1$ = ' num2str(L1(j))], 'FontSize', fnt);
    end

    % Std deviation
    yyaxis('right');
    hold on; j = 1;
    %plot(x(i, :), vs{j}.R(i, :)./vs{j}.R1(i, :) - 1, 'b--', 'LineWidth', 2);
    %plot(x(i, :), vs{j}.E(i, :)./vs{j}.R1(i, :) - 1, 'r--', 'LineWidth', 2);
    plot(x(i, :), sqrt(vs{j}.R1(i, :)), 'k.', 'LineWidth', 2, 'MarkerSize', 12);
    plot(x(i, :), sqrt(vs{j}.R(i, :)), 'b.', 'LineWidth', 2, 'MarkerSize', 12);
    plot(x(i, :), sqrt(vs{j}.E(i, :)), 'r.', 'LineWidth', 2, 'MarkerSize', 12);
    hold off; box off; grid off;
    h = gca; h.YColor = h.XColor;
    %ylabel(['$\sigma(X) | r = $ ' num2str(aj(i))], 'FontSize', fnt);
    if i == lena
        xlabel('$R_1 R_2^{-1}$', 'FontSize', fnt);
    end
end

for i = 1:lena
    % Plots of relative prob resurgence
    subplot(3, 2, 2*i);
    % Prob resurgence
    yyaxis('left');
    hold on; j = 2;
    %plot(x(i, :), Fs{j}.R(i, :)./Fs{j}.R1(i, :) - 1, 'b-', 'LineWidth', 2);
    %plot(x(i, :), Fs{j}.E(i, :)./Fs{j}.R1(i, :) - 1, 'r-', 'LineWidth', 2);
    plot(x(i, :), Fs{j}.R1(i, :), 'k-', 'LineWidth', 2);
    plot(x(i, :), Fs{j}.R(i, :), 'b-', 'LineWidth', 2);
    plot(x(i, :), Fs{j}.E(i, :), 'r-', 'LineWidth', 2);
    hold off; box off; grid off;
    h = gca; h.YColor = h.XColor;
    ymin = h.YLim(1); ylim([ymin 1.01]);
    %ylabel(['P$(X > 1) | r = $ ' num2str(aj(i))], 'FontSize', fnt);
    if i == 1
        title(['$\Lambda_1$ = ' num2str(L1(j))], 'FontSize', fnt);
    end

    % Std deviation
    yyaxis('right');
    hold on;
    %plot(x(i, :), vs{j}.R(i, :)./vs{j}.R1(i, :) - 1, 'b--', 'LineWidth', 2);
    %plot(x(i, :), vs{j}.E(i, :)./vs{j}.R1(i, :) - 1, 'r--', 'LineWidth', 2);
    plot(x(i, :), sqrt(vs{j}.R1(i, :)), 'k.', 'LineWidth', 2, 'MarkerSize', 12);
    plot(x(i, :), sqrt(vs{j}.R(i, :)), 'b.', 'LineWidth', 2, 'MarkerSize', 12);
    plot(x(i, :), sqrt(vs{j}.E(i, :)), 'r.', 'LineWidth', 2, 'MarkerSize', 12);
    hold off; box off; grid off;
    h = gca; h.YColor = h.XColor;
    ylabel(['$\sigma(X) | r = $ ' num2str(aj(i))], 'FontSize', fnt);
    if i == lena
        xlabel('$R_1 R_2^{-1}$', 'FontSize', fnt);
    end
end

