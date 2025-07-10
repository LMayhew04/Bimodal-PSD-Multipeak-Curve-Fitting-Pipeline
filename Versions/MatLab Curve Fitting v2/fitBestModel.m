function mdl = fitBestModel(tbl)
%FITBESTMODEL  Compare 2- and 3-component mixtures, choose lowest BIC.
%               • weighted least-squares (w = 1/x) boosts the fine tail
%               • multi-start (5 jittered guesses) avoids local minima
%
%   mdl = fitBestModel(tbl)   → struct with fields .type .BIC .pdf(x) .params

% ------------------------------------------------------------------------
x = tbl.ParticleDiameter;
y = tbl.SizeDistribution1;             % already area-normalised

xg = logspace(log10(min(x)),log10(max(x)),300).';
yg = interp1(x,y,xg,'linear','extrap');

cand = {@mNormLogn,@mNormRayl,@mNormNorm,@m3Comp}; % 3-component last
bestBIC = inf; best = struct;

for k = 1:numel(cand)
    c = cand{k}(xg,yg);
    if c.BIC < bestBIC
        bestBIC = c.BIC; best = c;
    end
end
mdl = best;
end

% ========================================================================
%                    2-COMPONENT CANDIDATES
% ========================================================================

function out = mNormLogn(x,y),  out = fitWith(@modelFcn,x,y,6,"Norm+Logn"); end
function out = mNormRayl(x,y),  out = fitWith(@modelFcn,x,y,5,"Norm+Rayl"); end
function out = mNormNorm(x,y),  out = fitWith(@modelFcn,x,y,6,"Norm+Norm"); end

function m = modelFcn(b,x,type) %#ok<INUSD>
switch type
    case "Norm+Logn"
        m = normalPDF(x,exp(b(1)),exp(b(2)),exp(b(3))) + ...
            lognPDF  (x,exp(b(4)),exp(b(5)),exp(b(6)));
    case "Norm+Rayl"
        m = normalPDF(x,exp(b(1)),exp(b(2)),exp(b(3))) + ...
            raylPDF  (x,exp(b(4)),exp(b(5)));
    case "Norm+Norm"
        m = normalPDF(x,exp(b(1)),exp(b(2)),exp(b(3))) + ...
            normalPDF(x,exp(b(4)),exp(b(5)),exp(b(6)));
end
end

% ========================================================================
%                    3-COMPONENT  (NEW)
% ========================================================================

function out = m3Comp(x,y)
    out = fitWith(@m3f,x,y,9,"Norm+Logn+Logn");
end
function m = m3f(b,x,~)
    m = normalPDF(x,exp(b(1)),exp(b(2)),exp(b(3))) + ...
        lognPDF  (x,exp(b(4)),exp(b(5)),exp(b(6))) + ...
        lognPDF  (x,exp(b(7)),exp(b(8)),exp(b(9)));
end

% ========================================================================
%                    CORE FITTER  (multi-start, weighted)
% ========================================================================

function out = fitWith(modelFun,x,y,p,type)
w   = 1./x;                       % emphasise fine tail
obj = @(b) rss(w.*modelFun(b,x,type), w.*y);

% ----- 5 jittered starts -------------------------------------------------
bestErr = inf; bestB = [];
b0 = log([median(x) std(x)/4 max(y) ...
          median(x)*6 0.6 max(y)   repmat([median(x)/10 0.4 max(y)*.1],1,max(0,p-6))]);

for s = 1:5
    guess = b0 .* exp(0.2*randn(size(b0)));
    [b,err] = fminsearch(obj,guess,optimset('Display','off','MaxIter',2e3));
    if err < bestErr, bestErr = err; bestB = b; end
end

yFit = modelFun(bestB,x,type);
yFit = yFit / trapz(x,yFit);     % ensure model area = 1

out.type   = type;
out.params = bestB;
out.pdf    = @(xq) modelFun(bestB,xq,type) / trapz(xq,modelFun(bestB,xq,type));
out.BIC    = bicScore(y,yFit,p);
end

% ========================================================================
%                         UTILITIES
% ========================================================================

function e = rss(a,b), e = sum((a-b).^2); end
function BIC = bicScore(y,yhat,k)
n=numel(y); rss=sum((y-yhat).^2);
BIC = n*log(rss/n)+k*log(n);
end