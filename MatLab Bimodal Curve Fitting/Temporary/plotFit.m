function plotFit(tbl, mdl, varargin)
%PLOTFIT  Plot PDF and CDF with histogram and fitted model.
%   plotFit(tbl, mdl, 'Title', ..., 'SaveDir', ...)

opt = struct('Title', '', 'SaveDir', 'Figures/');
opt = parseInputs(opt, varargin{:});

x = tbl.ParticleDiameter;
y = tbl.SizeDistribution1;

xg = mdl.x;
yfit = mdl.Fit;

% Axis limit (NEW): Max + 5% buffer, min at data min
xlimVals = [min(x), max(x) * 1.05];

figure('Visible', 'off'); clf;
subplot(2,1,1)
bar(x, y, 1, 'FaceAlpha', 0.5, 'DisplayName', 'Data');
hold on
plot(xg, yfit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit');
xlabel('Particle Diameter (μm)'); ylabel('Distribution (%)');
set(gca, 'XScale', 'log');
xlim(xlimVals);
legend; grid on;
if ~isempty(opt.Title), title(opt.Title); end

% CDF plot
subplot(2,1,2)
cdfData = cumsum(y) / sum(y);
cdfFit = cumsum(yfit) / sum(yfit);
plot(x, cdfData, 'bo-', 'DisplayName', 'Data CDF');
hold on
plot(xg, cdfFit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit CDF');
xlabel('Particle Diameter (μm)'); ylabel('Cumulative Fraction');
set(gca, 'XScale', 'log');
xlim(xlimVals);
legend; grid on;

if ~exist(opt.SaveDir, 'dir'), mkdir(opt.SaveDir); end
saveas(gcf, fullfile(opt.SaveDir, [opt.Title, '_fit.png']));
close(gcf);

end

function opt = parseInputs(opt, varargin)
for k = 1:2:numel(varargin)
    opt.(varargin{k}) = varargin{k+1};
end
end