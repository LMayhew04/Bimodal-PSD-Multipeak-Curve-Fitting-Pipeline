function plotFit(tbl, model, outPng)
%PLOTFIT  Create diagnostic plots (PDF, CDF, residuals) and save as PNG.
%         Adds residual stats and peak labels.

xFine  = logspace(log10(min(tbl.ParticleDiameter)), ...
                 log10(max(tbl.ParticleDiameter)), 400)';
pdfFit = model.pdf(xFine);

figure('Visible','off','Position',[100 100 800 600]);
tiledlayout(3,1,'TileSpacing','tight');

% --- PDF Overlay -----------------------------------------------------
nexttile;
semilogx(tbl.ParticleDiameter, tbl.SizeDistribution1, 'r.--', ...
         xFine, pdfFit, 'k-', 'LineWidth', 1.3);
ylabel('PDF'); legend('Data', 'Fit'); grid on;
title(sprintf('Peaks: %s', join(string(round(model.Peaks,2)), ', ')));

% --- CDF Overlay -----------------------------------------------------
nexttile;
cdfData = cumtrapz(tbl.ParticleDiameter, tbl.SizeDistribution1);
cdfFit  = cumtrapz(xFine, pdfFit);
cdfData = cdfData ./ cdfData(end);
cdfFit  = cdfFit  ./ cdfFit(end);
semilogx(tbl.ParticleDiameter, cdfData, 'r.--', ...
         xFine, cdfFit, 'k-'); ylabel('CDF'); grid on;

% --- Residuals -------------------------------------------------------
nexttile;
datInterp = interp1(tbl.ParticleDiameter, tbl.SizeDistribution1, ...
                    xFine, 'linear', 'extrap');
residual = datInterp - pdfFit;
semilogx(xFine, residual, 'b-'); grid on;
xlabel('Particle Diameter [\mum]');
ylabel('Residual');

% --- Title & Save ----------------------------------------------------
sgtitle(sprintf("%s | BIC = %.1f | MaxRes = %.3f | RMSE = %.3f", ...
    model.type, model.BIC, model.MaxResidual, model.RMSE));

[d,~,~] = fileparts(outPng);
if ~exist(d, 'dir'), mkdir(d); end
saveas(gcf, outPng); close(gcf);
end
