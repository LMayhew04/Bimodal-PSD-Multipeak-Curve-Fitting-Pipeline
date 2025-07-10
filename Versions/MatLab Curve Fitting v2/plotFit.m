function plotFit(tbl,model,outPng)
% Draw PDF, CDF, and residuals then save .png

    xFine  = logspace(log10(min(tbl.ParticleDiameter)),...
                      log10(max(tbl.ParticleDiameter)),400).';
    pdfFit = model.pdf(xFine);

    figure('Visible','off','Position',[100 100 800 600]);
    tiledlayout(3,1,'TileSpacing','tight')

    % --- PDF overlay -----------------------------------------------------
    nexttile
        semilogx(tbl.ParticleDiameter,tbl.SizeDistribution1,'r.--',...
                 xFine,pdfFit,'k-','LineWidth',1.3);
        ylabel('PDF'); legend('Data','Fit'); grid on

    % --- CDF overlay -----------------------------------------------------
    nexttile
        cdfData = cumtrapz(tbl.ParticleDiameter,tbl.SizeDistribution1);
        cdfFit  = cumtrapz(xFine,pdfFit);
        cdfData = cdfData ./ cdfData(end);
        cdfFit  = cdfFit  ./ cdfFit(end);
        semilogx(tbl.ParticleDiameter,cdfData,'r.--',xFine,cdfFit,'k-');
        ylabel('CDF'); grid on

    % --- Residuals -------------------------------------------------------
    nexttile
        datInterp = interp1(tbl.ParticleDiameter,tbl.SizeDistribution1,...
                            xFine,'linear','extrap');
        semilogx(xFine,datInterp-pdfFit,'b-'); grid on
        xlabel('Particle Diameter [\mum]'); ylabel('Residual')

    sgtitle(sprintf("%s  |  BIC = %.1f",model.type,model.BIC));

    % auto-create sub-folder and save
    [d,~,~] = fileparts(outPng);  if ~exist(d,'dir'), mkdir(d); end
    saveas(gcf,outPng);  close(gcf)
end