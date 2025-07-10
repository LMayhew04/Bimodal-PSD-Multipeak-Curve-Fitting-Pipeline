function runAllFits()
%RUNALLFITS  Batch-fit all Excel files in the Data/ folder.
%            Fits each sheet and saves PNGs and a summary CSV.

% --------- User Settings ------------------------------------------------
dataFolder = 'Data';
figureFolder = 'Figures';
customNames = {'ParticleDiameter', 'SizeDistribution1'}; % or your custom labels
weightExponent = 0.5;
sheetToUse = [2 4 6 8]; % Options: number, [2 3], 0, or 'all'

% --------- Setup Output File --------------------------------------------
summary = {};
summaryHeader = {'File', 'Sheet', 'Model', 'BIC', 'RMSE', 'MaxResidual'};

% --------- Loop Through Excel Files -------------------------------------
files = dir(fullfile(dataFolder, '*.xlsx'));
files = files(~startsWith({files.name}, '~$'));
files = files(~startsWith({files.name}, '~$'));  % Skip temp Excel lock files
for k = 1:length(files)
    [~, name, ~] = fileparts(files(k).name);
    filePath = fullfile(dataFolder, files(k).name);

    [~, sheets] = xlsfinfo(filePath);
    if isequal(sheetToUse, 'all')
        sheetList = 1:length(sheets);
    elseif isequal(sheetToUse, 0)
        sheetList = 2:length(sheets);
    elseif isnumeric(sheetToUse)
        sheetList = sheetToUse;
    else
        error('Invalid value for sheetToUse. Use numeric, 0, or ''all''.');
    end

    for s = sheetList
        try
            sheetName = sheets{s};  % extract actual sheet name
            tbl = loadParticleData(filePath, ...
                'Sheet', sheetName, ...
                'CustomNames', customNames);

            mdl = fitBestModel(tbl, ...
                'WeightExponent', weightExponent);

            figName = sprintf('%s_%s.png', name, sheetName);
            plotFit(tbl, mdl, fullfile(figureFolder, figName));

            summary(end+1, :) = {files(k).name, sheetName, ...
                             mdl.type, mdl.BIC, mdl.RMSE, mdl.MaxResidual}; %#ok<AGROW>
        catch ME
            warning('Failed: %s (Sheet: %s) — %s', files(k).name, sheets(s), ME.message);
        end
    end
end

% Write summary CSV if data exists
if ~isempty(summary)
    writetable(cell2table(summary, 'VariableNames', summaryHeader), ...
               'Fit_Summary.csv');
    disp('✅ All fittings completed successfully.');
    disp(['📝 Summary saved to: Fit_Summary.csv']);
    disp(['🖼️  Plots saved to folder: ' figureFolder]);
else
    warning('No summary results to write. All fits may have failed.');
end
end
