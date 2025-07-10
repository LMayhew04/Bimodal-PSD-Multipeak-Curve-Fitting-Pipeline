function runAllFits()
%RUNALLFITS  Batch-fit all Excel files in the Data/ folder.
%            Fits each sheet and saves PNGs and a summary CSV.
%
%   Update settings below to control sheet selection and fitting options.

% --------- User Settings ------------------------------------------------
dataFolder = 'Data';
figureFolder = 'Figures';
customNames = {'ParticleDiameter', 'SizeDistribution1'}; % or your custom labels
useThreeComponent = false;
weightExponent = 0.5;
sheetToUse = [2 3 4]; % set to 0 for all except 1st

% --------- Setup Output File --------------------------------------------
summary = {};
summaryHeader = {'File', 'Sheet', 'Model', 'BIC'};

% --------- Loop Through Excel Files -------------------------------------
files = dir(fullfile(dataFolder, '*.xls*'));
for k = 1:length(files)
    [~, name, ~] = fileparts(files(k).name);
    filePath = fullfile(dataFolder, files(k).name);

    [~, sheets] = xlsfinfo(filePath);
    if isequal(sheetToUse, 'all')
        sheetList = 1:length(sheets);           % all sheets
    elseif isequal(sheetToUse, 0)
        sheetList = 2:length(sheets);           % all except first
    elseif isnumeric(sheetToUse)
        sheetList = sheetToUse;                 % use specified sheets
    else
        error('Invalid value for sheetToUse. Use numeric, 0, or ''all''.');
    end

    for s = sheetList
        try
            tbl = loadParticleData(filePath, ...
                'Sheet', s, ...
                'CustomNames', customNames);

            mdl = fitBestModel(tbl, ...
                'WeightExponent', weightExponent);

            figName = sprintf('%s_%s.png', name, sheets{s});
            plotFit(tbl, mdl, fullfile(figureFolder, figName));

            summary(end+1, :) = {files(k).name, sheets{s}, mdl.type, mdl.BIC}; %#ok<AGROW>
        catch ME
            warning('Failed: %s (Sheet: %s) — %s', files(k).name, sheets{s}, ME.message);
        end
    end
end

% Write summary CSV
if ~isempty(summary)
    writetable(cell2table(summary, 'VariableNames', summaryHeader), ...
               'Fit_Summary.csv');
else
    warning('No summary results to write. All fits may have failed.');
end
end
