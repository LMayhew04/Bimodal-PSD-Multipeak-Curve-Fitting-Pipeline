function mdl = fitBestModel(tbl, varargin)
%FITBESTMODEL 2–4 component constrained fitting using fminsearch with penalties
% Works without Optimization Toolbox (no fmincon)

p = inputParser;
addParameter(p, 'WeightExponent', 0.5, @(a) isnumeric(a) && a >= 0);
parse(p, varargin{:});
opt = p.Results;

x = tbl.ParticleDiameter;
y = tbl.SizeDistribution1;
xg = logspace(-2, 2.5, 400)';  % fixed 0.01–300 µm
w = 1 ./ (xg .^ opt.WeightExponent);
w = w / max(w);
yg = interp1(x, y, xg, 'linear', 'extrap');

cand = {@m2Comp, @m3Comp, @m4Comp};
bestBIC = inf; best = struct;

for i = 1:numel(cand)
    c = cand{i}(xg, yg, w);
    if c.BIC < bestBIC
        bestBIC = c.BIC;
        best = c;
    end
end

% Diagnostics
best.MaxResidual = max(abs(best.Residual));
best.RMSE = sqrt(mean(best.Residual.^2));

% Bootstrapping
Nboot = 100;
boots = cell(Nboot,1);
for i = 1:Nboot
    idx = randsample(1:height(tbl), height(tbl), true);
    tboot = tbl(idx,:);
    try
        boots{i} = fitBestModel(tboot, 'WeightExponent', opt.WeightExponent);
    catch
        boots{i} = [];
    end
end
best.Bootstrap = boots(~cellfun('isempty',boots));
mdl = best;

% Models
function out = m2Comp(x, y, w)
    b0 = [2 0.4 0.03 10 0.6 0.03];
    fun = @(b) costWithBounds(b, x, y, w, @nlModel, ...
        [1 0.1 0 1 0.3 0], [100 2 1 100 2 1]);
    b = fminsearch(fun, b0, optimset('Display','off'));
    [yfit, peaks] = nlModel(b, x);
    out = wrap('Norm+Logn', b, yfit, y, x, w, peaks);
end

function out = m3Comp(x, y, w)
    b0 = [2 0.4 0.03 0.8 0.4 0.03 20 0.6 0.03];
    fun = @(b) costWithBounds(b, x, y, w, @nllModel, ...
        [1 0.1 0 0.2 0.2 0 5 0.3 0], [100 2 1 5 1 1 100 2 1]);
    b = fminsearch(fun, b0, optimset('Display','off'));
    [yfit, peaks] = nllModel(b, x);
    out = wrap('Norm+2×Logn', b, yfit, y, x, w, peaks);
end

function out = m4Comp(x, y, w)
    b0 = [2 0.4 0.03 0.8 0.4 0.03 20 0.6 0.03 0.4 0.2 0.01];
    fun = @(b) costWithBounds(b, x, y, w, @nlllModel, ...
        [1 0.1 0 0.2 0.2 0 5 0.3 0 0.3 0.1 0], ...
        [100 2 1 5 1 1 100 2 1 0.5 0.3 0.05]);
    b = fminsearch(fun, b0, optimset('Display','off'));
    [yfit, peaks] = nlllModel(b, x);
    out = wrap('Norm+3×Logn', b, yfit, y, x, w, peaks);
end

% Cost with penalty
function val = costWithBounds(b, x, y, w, f, lb, ub)
    penalty = 1e6 * sum(b < lb | b > ub);
    [yfit, ~] = f(b, x);
    val = sum(w .* (yfit - y).^2) + penalty;
end

function [yfit, peaks] = nlModel(b, x)
    p1 = normalPDF(x, b(1), b(2), b(3));
    p2 = lognPDF(x, b(4), b(5), b(6));
    yfit = p1 + p2; peaks = [b(1), b(4)];
end

function [yfit, peaks] = nllModel(b, x)
    p1 = normalPDF(x, b(1), b(2), b(3));
    p2 = lognPDF(x, b(4), b(5), b(6));
    p3 = lognPDF(x, b(7), b(8), b(9));
    yfit = p1 + p2 + p3; peaks = [b(1), b(4), b(7)];
end

function [yfit, peaks] = nlllModel(b, x)
    p1 = normalPDF(x, b(1), b(2), b(3));
    p2 = lognPDF(x, b(4), b(5), b(6));
    p3 = lognPDF(x, b(7), b(8), b(9));
    p4 = lognPDF(x, b(10), b(11), b(12));
    yfit = p1 + p2 + p3 + p4; peaks = [b(1), b(4), b(7), b(10)];
end

function out = wrap(type, b, yfit, ytrue, x, w, peaks)
    RSS = sum(w .* (yfit - ytrue).^2);
    n = numel(ytrue); k = numel(b);
    BIC = n*log(RSS/n) + k*log(n);
    out = struct('type', type, 'params', b, ...
        'pdf', @(xq) interp1(x, yfit, xq, 'linear', 'extrap'), ...
        'BIC', BIC, 'Residual', ytrue - interp1(x, yfit, x, 'linear', 'extrap'), ...
        'Peaks', peaks);
end
end
