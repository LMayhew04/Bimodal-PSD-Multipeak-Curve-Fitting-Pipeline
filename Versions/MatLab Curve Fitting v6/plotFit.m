function plotFit(tbl, model, outPng)
%PLOTFIT  Plot with y-axis clamping, spike filtering, and PDF smoothing

xFine = logspace(-2, 2.5, 400)';
pdfFit = model.pdf(xFine);

% Filter out spikes and clamp PDF
pdfFit(~isfinite(pdfFit)) = 0;  % remove NaN/Inf
pdfFit = max(pdfFit, 0);        % no negatives
pdfFit = min(pdfFit, 1);        % avoid extreme values
pdfFit = movmean(pdfFit, 3);    % smooth

pdfFit = pdfFit / trapz(xFine, pdfFit);  % normalize

% Bootstrap confidence bands
hasBoot = isfield(model, 'Bootstrap') && ~isempty(model.Bootstrap);
if hasBoot
    yMat = zeros(length(xFine), length(model.Bootstrap));
    for i = 1:length(model.Bootstrap)
        y = model.Bootstrap{i}.pdf(xFine);
        y(~isfinite(y)) = 0;
        y = max(min(y, 1), 0);
        yMat(:, i) = movmean(y, 3) / trapz(xFine, movmean(y, 3));
    end
    lower = prctile(yMat, 2.5, 2);
    upper = prctile(yMat, 97.5, 2);
end

figure('Visible','off','Position',[100 100 800 600]);
tiledlayout(3,1,'TileSpacing','tight');

% --- PDF -----------------------------------------------------
nexttile; hold on;
if hasBoot
    fill([xFine; flipud(xFine)], [lower; flipud(upper)], ...
        [0.9 0.9 1], 'EdgeColor','none', 'FaceAlpha', 0.5);
end
semilogx(tbl.ParticleDiameter, tbl.SizeDistribution1, 'r.--');
semilogx(xFine, pdfFit, 'k-', 'LineWidth', 1.3);
for i = 1:length(model.Peaks)
    xpk = model.Peaks(i);
    [~, idx] = min(abs(xFine - xpk));
    ypk = pdfFit(idx);
    if isfinite(ypk)
        text(xpk, ypk, sprintf('%.2f', ypk), ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','center');
    end
end
xmax = max(tbl.ParticleDiameter(tbl.SizeDistribution1 > 0)) * 1.1;
xlim([min(xFine), xmax]); ylim([0, 0.1]);
ylabel('PDF'); legend('Data', 'Fit'); grid on; hold off;
title(sprintf('Peaks: %s', join(string(round(model.Peaks,2)), ', ')));

% --- CDF -----------------------------------------------------
nexttile;
cdfData = cumtrapz(tbl.ParticleDiameter, tbl.SizeDistribution1);
cdfFit = cumtrapz(xFine, pdfFit);
cdfData = cdfData / max(cdfData);
cdfFit  = cdfFit  / max(cdfFit);
semilogx(tbl.ParticleDiameter, cdfData, 'r.--', xFine, cdfFit, 'k-');
ylabel('CDF'); grid on;

% --- Residual ------------------------------------------------
nexttile;
datInterp = interp1(tbl.ParticleDiameter, tbl.SizeDistribution1, ...
                    xFine, 'linear', 'extrap');
residual = datInterp - pdfFit;
semilogx(xFine, residual, 'b-'); hold on; grid on;
xlabel('Particle Diameter [\mum]'); ylabel('Residual');
thresh = 0.005;
outliers = abs(residual) > thresh;
plot(xFine(outliers), residual(outliers), 'ro', 'MarkerSize', 4);

sgtitle(sprintf("%s | BIC = %.1f | MaxRes = %.3f | RMSE = %.3f", ...
    model.type, model.BIC, model.MaxResidual, model.RMSE));

[d,~,~] = fileparts(outPng);
if ~exist(d, 'dir'), mkdir(d); end
saveas(gcf, outPng); close(gcf);
end
