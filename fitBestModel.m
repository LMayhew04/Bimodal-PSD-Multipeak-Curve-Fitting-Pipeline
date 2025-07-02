function model = fitBestModel(diam, freq, opts)
% Fit a mixture of lognormal and/or normal components with explicit
% penalty for missed data peaks and minimum width constraint.

if nargin < 3
    opts = struct;
end
opts = setdefault(opts, 'numComp', 3);
opts = setdefault(opts, 'type', 'mix'); % 'logn', 'norm', or 'mix'
opts = setdefault(opts, 'minWidth', 0.03*range(diam)); % or by log
opts = setdefault(opts, 'maxTries', 20);
opts = setdefault(opts, 'peakTol', 0.07*range(diam));

bestCost = Inf; bestParams = [];
N = opts.numComp;

for attempt = 1:opts.maxTries
    mus = min(diam) + (max(diam)-min(diam)) * rand(1,N);
    sigs = repmat(0.15*range(diam),1,N) .* abs(1+0.1*randn(1,N));
    ws = rand(1,N); ws = ws/sum(ws);
    params0 = reshape([mus; sigs; ws],1,[]);
    
    costfun = @(p) pdfCostWithPenalties(p, diam, freq, opts);

    options = optimset('Display','off', 'MaxIter', 1200, 'TolFun',1e-7, 'TolX',1e-6);
    [pfit, fval] = fminsearch(costfun, params0, options);
    if fval < bestCost
        bestCost = fval;
        bestParams = pfit;
    end
end

model = buildModelFromParams(bestParams, opts, diam);
end

function cost = pdfCostWithPenalties(params, diam, freq, opts)
    N = opts.numComp;
    mus = params(1:3:end);
    sigs = abs(params(2:3:end));
    ws = abs(params(3:3:end)); ws = ws/sum(ws);

    xq = linspace(min(diam), max(diam), 400);
    y = zeros(size(xq));
    for i = 1:N
        if strcmp(opts.type, "logn") || strcmp(opts.type, 'logn')
            y = y + ws(i)*lognPDF(xq, mus(i), sigs(i));
        elseif strcmp(opts.type, "norm") || strcmp(opts.type, 'norm')
            y = y + ws(i)*normalPDF(xq, mus(i), sigs(i));
        else % 'mix': alternate logn/norm by component
            if mod(i,2)==1
                y = y + ws(i)*lognPDF(xq, mus(i), sigs(i));
            else
                y = y + ws(i)*normalPDF(xq, mus(i), sigs(i));
            end
        end
    end

    if any(~isfinite(xq)) || any(~isfinite(y)) || any(~isfinite(diam))
        cost = 1e12;
        return
    end
    y_at_data = interp1(xq, y, diam, 'linear', 'extrap');
    if any(~isfinite(y_at_data))
        cost = 1e12;
        return
    end
    y_at_data(y_at_data<0) = 0;
    mse = mean((y_at_data - freq).^2);

    dataLocs = simplePeaks(diam, freq, max(freq)*0.05);
    modelLocs = simplePeaks(xq, y, max(y)*0.04);

    missed = 0;
    for d = dataLocs(:)'
        if all(abs(modelLocs - d) > opts.peakTol)
            missed = missed + 1;
        end
    end
    peakPenalty = missed * 1e3;
    widthPenalty = sum(sigs < opts.minWidth) * 2e3;

    cost = mse + peakPenalty + widthPenalty;
end

function model = buildModelFromParams(params, opts, diam)
    N = opts.numComp;
    mus = params(1:3:end);
    sigs = abs(params(2:3:end));
    ws = abs(params(3:3:end)); ws = ws/sum(ws);
    model.pdf = @(x) mixturePDF(x, mus, sigs, ws, opts.type);
    model.params = params;
    model.components = struct('mu',mus,'sigma',sigs,'weight',ws);
    model.type = opts.type;
    model.domain = [min(diam), max(diam)];
end

function y = mixturePDF(x, mus, sigs, ws, type)
    N = numel(mus);
    y = zeros(size(x));
    for i = 1:N
        if strcmp(type, "logn") || strcmp(type, 'logn')
            y = y + ws(i)*lognPDF(x, mus(i), sigs(i));
        elseif strcmp(type, "norm") || strcmp(type, 'norm')
            y = y + ws(i)*normalPDF(x, mus(i), sigs(i));
        else
            if mod(i,2)==1
                y = y + ws(i)*lognPDF(x, mus(i), sigs(i));
            else
                y = y + ws(i)*normalPDF(x, mus(i), sigs(i));
            end
        end
    end
end

function y = lognPDF(x, mu, sigma)
    y = zeros(size(x));
    positive = x > 0;
    y(positive) = (1./(x(positive)*sigma*sqrt(2*pi))) .* exp(-((log(x(positive))-log(mu)).^2)/(2*sigma^2));
end

function y = normalPDF(x, mu, sigma)
    y = (1/(sigma*sqrt(2*pi))) * exp(-(x-mu).^2/(2*sigma^2));
end

function s = setdefault(s, field, val)
    if ~isfield(s, field)
        s.(field) = val;
    end
end

function peakLocs = simplePeaks(x, y, minProm)
if nargin < 3, minProm = 0; end
peakLocs = [];
for k = 2:length(y)-1
    if y(k) > y(k-1) && y(k) > y(k+1) && y(k) > minProm
        peakLocs(end+1) = x(k); %#ok<AGROW>
    end
end
end