function mdl = fitBestModel(tbl, varargin)
%FITBESTMODEL  Automatically select best 2- or 3-component model based on BIC.
%   mdl = fitBestModel(tbl, 'WeightExponent', 0.5)
%
%   Tries:
%     - Normal + Lognormal
%     - Normal + Normal
%     - Normal + 2×Lognormal
%
%   Selects the model with the lowest BIC and computes residual stats.

% Parse options
p = inputParser;
addParameter(p, 'WeightExponent', 0.5, @(a) isnumeric(a) && a >= 0);
parse(p, varargin{:});
opt = p.Results;

x = tbl.ParticleDiameter;
y = tbl.SizeDistribution1;

% Log-spaced interpolation grid
xg = logspace(log10(min(x)), log10(max(x)), 300)';
yg = interp1(x, y, xg, 'linear', 'extrap');

% Weight function
w = 1 ./ (xg .^ opt.WeightExponent);
w = w / max(w);  % Normalize

% Candidates
cand = {@mNormLogn, @mNormNorm, @m3Comp};
bestBIC = inf; best = struct;

for k = 1:numel(cand)
    c = cand{k}(xg, yg, w);
    if c.BIC < bestBIC
        bestBIC = c.BIC;
        best = c;
    end
end

% Attach residual stats
best.MaxResidual = max(abs(best.Residual));
best.RMSE = sqrt(mean(best.Residual .^ 2));
mdl = best;

% ---------------------- Candidate Models -------------------------------

function out = mNormLogn(x, y, w)
    beta0 = [log(median(x)) log(std(x)/4) log(max(y)) ...
             log(median(x)*6) log(0.6) log(max(y))];
    obj = @(b) weightedRSS(b, x, y, w, @nlModel);
    bEst = fminsearch(obj, beta0, optimset('Display','off'));
    [yFit, pk] = nlModel(bEst, x);
    out = package('Norm+Logn', bEst, yFit, y, x, w, pk);
end

function out = mNormNorm(x, y, w)
    beta0 = [log(median(x)) log(std(x)/4) log(max(y)) ...
             log(median(x)*2) log(0.4) log(max(y))];
    obj = @(b) weightedRSS(b, x, y, w, @nnModel);
    bEst = fminsearch(obj, beta0, optimset('Display','off'));
    [yFit, pk] = nnModel(bEst, x);
    out = package('Norm+Norm', bEst, yFit, y, x, w, pk);
end

function out = m3Comp(x, y, w)
    beta0 = [log(median(x)) log(std(x)/4) log(max(y)/2) ...
             log(median(x)/2) log(0.4) log(max(y)/3) ...
             log(median(x)*2) log(0.6) log(max(y)/3)];
    obj = @(b) weightedRSS(b, x, y, w, @nllModel);
    bEst = fminsearch(obj, beta0, optimset('Display','off'));
    [yFit, pk] = nllModel(bEst, x);
    out = package('Norm+2×Logn', bEst, yFit, y, x, w, pk);
end

% ------------------- Model Equations ------------------------

function [yFit, peaks] = nlModel(b, x)
    p1 = normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3)));
    p2 = lognPDF(x, exp(b(4)), exp(b(5)), exp(b(6)));
    yFit = p1 + p2;
    peaks = [exp(b(1)), exp(b(4))];
end

function [yFit, peaks] = nnModel(b, x)
    p1 = normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3)));
    p2 = normalPDF(x, exp(b(4)), exp(b(5)), exp(b(6)));
    yFit = p1 + p2;
    peaks = [exp(b(1)), exp(b(4))];
end

function [yFit, peaks] = nllModel(b, x)
    p1 = normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3)));
    p2 = lognPDF(x, exp(b(4)), exp(b(5)), exp(b(6)));
    p3 = lognPDF(x, exp(b(7)), exp(b(8)), exp(b(9)));
    yFit = p1 + p2 + p3;
    peaks = [exp(b(1)), exp(b(4)), exp(b(7))];
end

% ----------------- Utilities -------------------------

function val = weightedRSS(b, x, y, w, modelFunc)
    [yFit, ~] = modelFunc(b, x);
    val = sum(w .* (yFit - y).^2);
end

function out = package(type, bEst, yFit, yTrue, x, w, peaks)
    RSS = sum(w .* (yFit - yTrue).^2);
    n = numel(yTrue); k = numel(bEst);
    BIC = n*log(RSS/n) + k*log(n);
    out = struct('type', type, 'params', bEst, ...
                 'pdf', @(xq) interp1(x, yFit, xq, 'linear', 'extrap'), ...
                 'BIC', BIC, ...
                 'Residual', yTrue - interp1(x, yFit, x, 'linear', 'extrap'), ...
                 'Peaks', peaks);
end
end
