% RUNALLFITS - Run full pipeline on all data sheets and files
% Centralizes all user settings.

% SETTINGS
sheetMode = 'allExceptFirst'; % Options: 'all', 'allExceptFirst', {'SheetName1', ...}, or 'single'
dataDir = 'Data/';
figDir = 'Figures/';
summaryFile = 'Fit_Summary.csv';

if exist(summaryFile, 'file'), delete(summaryFile); end

dataFiles = dir(fullfile(dataDir, '*.xlsx'));
dataFiles = dataFiles(~startsWith({dataFiles.name}, '~$')); % Skip temp files

summary = {};

for f = 1:length(dataFiles)
    [sheets, ~] = xlsfinfo(fullfile(dataDir, dataFiles(f).name));
    if strcmp(sheetMode, 'allExceptFirst')
        sheetList = sheets(2:end);
    elseif iscell(sheetMode)
        sheetList = sheetMode;
    elseif strcmp(sheetMode, 'single')
        sheetList = sheets(1);
    else
        sheetList = sheets;
    end
    sheetList = cellstr(sheetList); % <-- Ensures always cell array (key fix)
    for s = 1:length(sheetList)
        try
            tbl = loadParticleData(fullfile(dataDir, dataFiles(f).name), sheetList{s});
            mdl = fitBestModel(tbl);
            plotFit(tbl, mdl, ...
                'Title', sprintf('%s - %s', dataFiles(f).name, sheetList{s}), ...
                'SaveDir', figDir);
            summary(end+1, :) = {dataFiles(f).name, sheetList{s}, ...
                mdl.Model, mdl.BIC, mdl.RMSE, mdl.MaxResidual, ...
                numel(mdl.PeakLoc)};
        catch err
            warning('Failed on %s:%s - %s', dataFiles(f).name, sheetList{s}, err.message);
        end
    end
end

% Write summary CSV
header = {'File', 'Sheet', 'Model', 'BIC', 'RMSE', 'MaxResidual', 'nPeaks'};
writecell([header; summary], summaryFile);