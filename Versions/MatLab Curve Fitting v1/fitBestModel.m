function mdl = fitBestModel(tbl)
%FITBESTMODEL  Pick best 2-peak model (Norm+Logn | Norm+Rayl | Norm+Norm)
%              using fminsearch optimisation.
% ------------------------------------------------------------------------

xdata = tbl.ParticleDiameter;
ydata = tbl.SizeDistribution1;

xgrid = logspace(log10(min(xdata)),log10(max(xdata)),300).';
ygrid = interp1(xdata,ydata,xgrid,'linear','extrap');

% run three candidates
fn   = {@modelNormLogn,@modelNormRayl,@modelNormNorm};
best = []; bestBIC = inf;

for k = 1:numel(fn)
    cand = fn{k}(xgrid,ygrid);
    if cand.BIC < bestBIC
        bestBIC = cand.BIC;
        best    = cand;
    end
end
mdl = best;
end   % ================= end MAIN function ===============================

%%                       CANDIDATE MODEL FITS
%   We optimise log-sigma & log-gain so parameters stay positive.

function out = modelNormLogn(x,y)
% Normal + Lognormal  (6 params)
beta0 = [log(median(x))  log(std(x)/4)  log(max(y))  ...
         log(median(x)*6) log(0.6)       log(max(y))];

obj = @(b) rss( ...
       normalPDF(x,exp(b(1)),exp(b(2)),exp(b(3))) + ...
       lognPDF  (x,exp(b(4)),exp(b(5)),exp(b(6))), y);

bEst = fminsearch(obj,beta0,optimset('Display','off'));

yFit = normalPDF(x,exp(bEst(1)),exp(bEst(2)),exp(bEst(3))) + ...
       lognPDF  (x,exp(bEst(4)),exp(bEst(5)),exp(bEst(6)));

out.type   = "Norm+Logn";
out.params = bEst;
out.pdf    = @(xq) normalPDF(xq,exp(bEst(1)),exp(bEst(2)),exp(bEst(3))) + ...
                   lognPDF  (xq,exp(bEst(4)),exp(bEst(5)),exp(bEst(6)));
out.BIC    = bicScore(y,yFit,numel(bEst));
end
% ------------------------------------------------------------------------
function out = modelNormRayl(x,y)
beta0 = [log(median(x))  log(std(x)/4)  log(max(y)) ...
         log(median(x)*3)               log(max(y))];

obj = @(b) rss( ...
       normalPDF(x,exp(b(1)),exp(b(2)),exp(b(3))) + ...
       raylPDF  (x,exp(b(4)),exp(b(5))), y);

bEst = fminsearch(obj,beta0,optimset('Display','off'));

yFit = normalPDF(x,exp(bEst(1)),exp(bEst(2)),exp(bEst(3))) + ...
       raylPDF  (x,exp(bEst(4)),exp(bEst(5)));

out.type = "Norm+Rayl"; out.params = bEst;
out.pdf  = @(xq) normalPDF(xq,exp(bEst(1)),exp(bEst(2)),exp(bEst(3))) + ...
                  raylPDF (xq,exp(bEst(4)),exp(bEst(5)));
out.BIC  = bicScore(y,yFit,numel(bEst));
end
% ------------------------------------------------------------------------
function out = modelNormNorm(x,y)
beta0 = [log(median(x)/3) log(std(x)/5)  log(max(y)/2) ...
         log(median(x)*3) log(std(x))    log(max(y)/2)];

obj = @(b) rss( ...
       normalPDF(x,exp(b(1)),exp(b(2)),exp(b(3))) + ...
       normalPDF(x,exp(b(4)),exp(b(5)),exp(b(6))), y);

bEst = fminsearch(obj,beta0,optimset('Display','off'));

yFit = normalPDF(x,exp(bEst(1)),exp(bEst(2)),exp(bEst(3))) + ...
       normalPDF(x,exp(bEst(4)),exp(bEst(5)),exp(bEst(6)));

out.type = "Norm+Norm"; out.params = bEst;
out.pdf  = @(xq) normalPDF(xq,exp(bEst(1)),exp(bEst(2)),exp(bEst(3))) + ...
                  normalPDF(xq,exp(bEst(4)),exp(bEst(5)),exp(bEst(6)));
out.BIC  = bicScore(y,yFit,numel(bEst));
end

% ========================================================================
%                             Helpers
% ========================================================================

function e = rss(model, data)      % residual-sum-of-squares
e = sum((model - data).^2);
end
function BIC = bicScore(y,yhat,k)
n = numel(y); rss = sum((y-yhat).^2);
BIC = n*log(rss/n) + k*log(n);
end