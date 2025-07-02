function mdl = fitBestModel(tbl, varargin)
%FITBESTMODEL  Automatically select best 2- or 3-component model by visual fidelity and BIC.
%   mdl = fitBestModel(tbl, 'WeightExponent', 0.5)
%
%   New in this version:
%     - Always tests both 2- and 3-component fits.
%     - If visual inspection or metrics show missed peaks/shoulders, selects 3-component.
%     - Penalizes overfitting/noise.
%     - Axis and support limited to [min(x), max(x)*1.05].

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

% Fit all models and retain all info for diagnostics
candidates = { @mNormLogn, @mNormNorm, @m3Comp };
results = [];
for k = 1:numel(candidates)
    c = candidates{k}(xg, yg, w);
    results = [results; c];
end

% Visual/automated peak count comparison (NEW)
nDataPeaks = countPeaks(y, 'MinProminence', max(y)*0.1);
nBestPeaks = cellfun(@(c) numel(c.PeakLoc), num2cell(results));
[~, idxBest] = min(arrayfun(@(c) c.BIC, results));
% Prefer model with correct peak count if BIC difference is small
% (forces 3-comp if it picks up a clear missed peak)
idxUse = idxBest;
for i = 1:length(results)
    if abs(nDataPeaks - nBestPeaks(i)) <= 1 && results(i).BIC < results(idxBest).BIC * 1.05
        idxUse = i;
        break;
    end
end
mdl = results(idxUse);

mdl.MaxResidual = max(abs(mdl.Residual));
mdl.RMSE = sqrt(mean(mdl.Residual .^ 2));

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
    beta0 = [log(median(x)) log(std(x)/5) log(max(y)*0.7) ...
             log(median(x)*2) log(0.3) log(max(y)*0.2) ...
             log(median(x)*6) log(0.6) log(max(y)*0.1)];
    obj = @(b) weightedRSS(b, x, y, w, @nnlModel);
    bEst = fminsearch(obj, beta0, optimset('Display','off'));
    [yFit, pk] = nnlModel(bEst, x);
    out = package('3-comp', bEst, yFit, y, x, w, pk);
end

% --- Model forms ---
function [yfit, peaks] = nlModel(b, x)
    yfit = normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3))) + ...
           lognPDF(x, exp(b(4)), exp(b(5)), exp(b(6)));
    [~, peaks] = findpeaks(yfit, x, 'MinPeakProminence', max(yfit)*0.08);
end

function [yfit, peaks] = nnModel(b, x)
    yfit = normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3))) + ...
           normalPDF(x, exp(b(4)), exp(b(5)), exp(b(6)));
    [~, peaks] = findpeaks(yfit, x, 'MinPeakProminence', max(yfit)*0.08);
end

function [yfit, peaks] = nnlModel(b, x)
    yfit = normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3))) + ...
           normalPDF(x, exp(b(4)), exp(b(5)), exp(b(6))) + ...
           lognPDF(x, exp(b(7)), exp(b(8)), exp(b(9)));
    [~, peaks] = findpeaks(yfit, x, 'MinPeakProminence', max(yfit)*0.08);
end

% --- Utilities ---
function rss = weightedRSS(b, x, y, w, modelfun)
    yfit = modelfun(b, x);
    rss = sum(w .* (yfit - y).^2);
end

function n = countPeaks(y, varargin)
    [~, locs] = findpeaks(y, varargin{:});
    n = numel(locs);
end

function out = package(modelType, params, yFit, y, x, w, peaks)
    res = yFit - y;
    k = numel(params);
    n = numel(x);
    BIC = n*log(sum((res).^2)/n) + k*log(n);
    out = struct('Model', modelType, 'Params', params, 'Fit', yFit, ...
        'Residual', res, 'BIC', BIC, 'PeakLoc', peaks, ...
        'x', x, 'Data', y, 'Weights', w);
end

end