function mdl = fitBestModel(tbl, varargin)
%FITBESTMODEL  Fit 2- or 3-component mixture model to particle data.
%   mdl = fitBestModel(tbl) fits models to the input table with columns:
%       - ParticleDiameter
%       - SizeDistribution1 (area-normalized PDF)
%
%   Optional name-value pairs:
%       'UseThreeComponent' - true/false (default: false)
%       'WeightExponent'    - exponent alpha in w = 1 / x^alpha (default: 0.5)
%
%   Output 'mdl' contains fields: .type, .BIC, .pdf, .params

% Parse options
p = inputParser;
addParameter(p, 'UseThreeComponent', false, @islogical);
addParameter(p, 'WeightExponent', 0.5, @(a) isnumeric(a) && a >= 0);
parse(p, varargin{:});
opt = p.Results;

x = tbl.ParticleDiameter;
y = tbl.SizeDistribution1;

% Interpolate to common log-spaced grid
xg = logspace(log10(min(x)), log10(max(x)), 300)';
yg = interp1(x, y, xg, 'linear', 'extrap');

% Apply weighting: w = 1 / x^alpha
w = 1 ./ (xg .^ opt.WeightExponent);
w = w / max(w);  % Normalize weights

% Select candidate models
cand = {@mNormLogn, @mNormNorm};  % Always include these two
if opt.UseThreeComponent
    cand{end+1} = @m3Comp;        % Optional 3-component model
end

% Try each model and keep the best by BIC
bestBIC = inf; best = struct;
for k = 1:numel(cand)
    c = cand{k}(xg, yg, w);
    if c.BIC < bestBIC
        bestBIC = c.BIC;
        best = c;
    end
end

mdl = best;

% ======================= Model Definitions ==========================

function out = mNormLogn(x, y, w)
    % Normal + Lognormal
    beta0 = [log(median(x)) log(std(x)/4) log(max(y)) ...
             log(median(x)*6) log(0.6) log(max(y))];
    obj = @(b) sum(w .* ( ...
        normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3))) + ...
        lognPDF(x, exp(b(4)), exp(b(5)), exp(b(6))) - y ).^2);
    bEst = fminsearch(obj, beta0, optimset('Display','off'));
    yFit = normalPDF(x, exp(bEst(1)), exp(bEst(2)), exp(bEst(3))) + ...
           lognPDF(x, exp(bEst(4)), exp(bEst(5)), exp(bEst(6)));
    out = makeOutput('Norm+Logn', bEst, yFit, y, x, w);
end

function out = mNormNorm(x, y, w)
    % Normal + Normal
    beta0 = [log(median(x)) log(std(x)/4) log(max(y)) ...
             log(median(x)*2) log(0.4) log(max(y))];
    obj = @(b) sum(w .* ( ...
        normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3))) + ...
        normalPDF(x, exp(b(4)), exp(b(5)), exp(b(6))) - y ).^2);
    bEst = fminsearch(obj, beta0, optimset('Display','off'));
    yFit = normalPDF(x, exp(bEst(1)), exp(bEst(2)), exp(bEst(3))) + ...
           normalPDF(x, exp(bEst(4)), exp(bEst(5)), exp(bEst(6)));
    out = makeOutput('Norm+Norm', bEst, yFit, y, x, w);
end

function out = m3Comp(x, y, w)
    % Normal + Lognormal + Lognormal
    beta0 = [log(median(x)) log(std(x)/4) log(max(y)/2) ...
             log(median(x)/2) log(0.4) log(max(y)/3) ...
             log(median(x)*2) log(0.6) log(max(y)/3)];
    obj = @(b) sum(w .* ( ...
        normalPDF(x, exp(b(1)), exp(b(2)), exp(b(3))) + ...
        lognPDF(x, exp(b(4)), exp(b(5)), exp(b(6))) + ...
        lognPDF(x, exp(b(7)), exp(b(8)), exp(b(9))) - y ).^2);
    bEst = fminsearch(obj, beta0, optimset('Display','off'));
    yFit = normalPDF(x, exp(bEst(1)), exp(bEst(2)), exp(bEst(3))) + ...
           lognPDF(x, exp(bEst(4)), exp(bEst(5)), exp(bEst(6))) + ...
           lognPDF(x, exp(bEst(7)), exp(bEst(8)), exp(bEst(9)));
    out = makeOutput('Norm+2×Logn', bEst, yFit, y, x, w);
end

function out = makeOutput(type, bEst, yFit, yTrue, x, w)
    RSS = sum(w .* (yFit - yTrue).^2);
    n = numel(yTrue);
    k = numel(bEst);
    BIC = n*log(RSS/n) + k*log(n);
    out = struct('type', type, 'params', bEst, ...
                 'pdf', @(xq) interp1(x, yFit, xq, 'linear', 'extrap'), ...
                 'BIC', BIC);
end
end
