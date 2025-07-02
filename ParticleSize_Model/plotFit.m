function plotFit(diam, pdf_data, model, outPng)
% plotFit - v4-style overlays, semilogx, CDF, residuals, tight axis scaling

if isempty(diam) || all(~isfinite(diam))
    warning('No finite diameter data found! Plot will be skipped.');
    return
end

min_d = min(diam);
max_d = max(diam);
range_d = max_d - min_d;
buffer = 0.05 * range_d;
if buffer < 1e-5, buffer = 1; end
x_lo = max(0, min_d - buffer);
x_hi = max_d + buffer;

x_plot = linspace(x_lo, x_hi, 400);
y_model = model.pdf(x_plot);

figure('Visible','off');

% --- Subplot 1: Empirical PDF and Model Overlay (semilogx) ---
subplot(3,1,1);
semilogx(diam, pdf_data, 'r--', 'LineWidth', 2); hold on;
semilogx(x_plot, y_model, 'k-', 'LineWidth', 2);
xlabel('Particle Diameter (\mum)');
ylabel('PDF');
title('Particle Size Distribution & Model Fit');
legend('Data (PDF)', 'Model', 'Location', 'Best');
grid on;
limits = sort([x_lo, x_hi]);
xlim(limits);

% --- Subplot 2: CDF Overlay ---
subplot(3,1,2);
cdf_data = cumsum(pdf_data) / sum(pdf_data);
cdf_model = cumsum(y_model) / sum(y_model);
semilogx(diam, cdf_data, 'ro--', 'LineWidth', 1.4, 'MarkerSize', 4); hold on;
semilogx(x_plot, cdf_model, 'k-', 'LineWidth', 2);
ylabel('CDF');
legend('Data (CDF)', 'Model', 'Location', 'Best');
grid on;
xlim(limits);

% --- Subplot 3: Residuals ---
subplot(3,1,3);
model_at_data = interp1(x_plot, y_model, diam, 'linear', 'extrap');
residuals = pdf_data - model_at_data;
semilogx(diam, residuals, 'bo-', 'LineWidth', 1.2, 'MarkerFaceColor', 'b');
xlabel('Particle Diameter (\mum)');
ylabel('Residual (Data - Model)');
title('Fit Residuals');
grid on;
xlim(limits);

rmse = sqrt(mean(residuals.^2));
max_resid = max(abs(residuals));
mae = mean(abs(residuals));
fit_percent_error = 100 * mae / max(pdf_data);

text_x = limits(1) + 0.65*(limits(2) - limits(1));
text_y = max(residuals)*0.85;
txt = sprintf('RMSE: %.3g\nMax resid: %.3g\nMean abs err: %.3g (%.2f%%)', ...
    rmse, max_resid, mae, fit_percent_error);
subplot(3,1,3);
text(text_x, text_y, txt, 'FontSize', 10, 'Color', [0.2 0.2 0.2]);

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, outPng);
close(gcf);
end