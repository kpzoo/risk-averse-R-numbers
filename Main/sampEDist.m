% Sample from distributions of Rj to compute E
function [RmE, RlE, RhE, p1E, Esamp] = sampEDist(nsamps, nday, distR, nDeme, istarts, Rgrid)

% E optimal distribution by sampling 
RmE = zeros(1, nday); RlE = RmE; RhE = RmE; 

% At every time sample from Rj distributions and compute E
for i = 1:nday
    
    % Individual distribution samples
    xDeme = zeros(nDeme, nsamps);
    for j = 1:nDeme
        if i >= istarts(j)
            xDeme(j, :) = datasample(Rgrid, nsamps, 'Weights', distR{j}(i-istarts(j)+1, :));
        end
    end

    % E optimal samples for this day
    Esamp = mean(xDeme.^2)./mean(xDeme);
    % Prob of E > 1
    p1E(i) = length(Esamp(Esamp > 1))/nsamps;

    % Statistics of E designs
    RmE(i) = mean(Esamp);
    Equants = quantile(Esamp, [0.025, 0.975]);
    RlE(i) = Equants(1); RhE(i) = Equants(2);
end
