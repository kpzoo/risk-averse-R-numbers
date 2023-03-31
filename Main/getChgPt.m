% Find changepoints for reproduction numbers and prob > 1
function [tR, t1] = getChgPt(t, R, p1)

% Assumptions and notes
% - change times for transition of R below to above 1
% - change times for P(R>1) below to above 0.95

% Normalise about thresholds
R = R - 1; p1 = p1 - 0.95;

% Find times for reproduction number to transition
tR = t(R(2:end) > 0 & R(1:end-1) < 0);
t1 = t(p1(2:end) > 0 & p1(1:end-1) < 0);

