%% runAllFits.m  –  Batch-fit particle-size data
%   ▸ SECTION A  : process **Sheet 2 only** for every workbook
%   ▸ SECTION B  : process **all sheets EXCEPT Sheet 1**
% Toggle by setting the corresponding if-flag to true and the other to false.

clear; clc;

rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir,'Data');
figDir  = fullfile(rootDir,'Figures');

addpath(fullfile(rootDir,'pdfHelpers')); % make sure PDF helpers are seen

files   = dir(fullfile(dataDir,'**','*.xls*')); % recurse sub-folders
results = {};    % collect rows as a cell array

%%%  SECTION A: SHEET 2 ONLY %%%

if false           % set to false when you want Section B instead
    wb = waitbar(0,'Processing files (Sheet 2 only)…');
    for i = 1:numel(files)
        waitbar(i/numel(files),wb,files(i).name)
        try
            fullName = fullfile(files(i).folder,files(i).name);

            % Sheet argument omitted – loadParticleData default is Sheet 2
            tbl = loadParticleData(fullName);

            mdl = fitBestModel(tbl);

            base   = erase(files(i).name,{'.xlsx','.xls'});
            pngName = fullfile(figDir,[base '.png']);
            plotFit(tbl,mdl,pngName);

            results = [results;
                       {files(i).name + ":Sheet2", mdl.type, mdl.BIC}];
        catch ME
            warning("Skipped %s (Sheet 2): %s",files(i).name,ME.message);
        end
    end
    close(wb)
end

%%% SECTION B: EVERY SHEET EXCEPT THE FIRST %%%

if true          % set to true when you want multi-sheet processing
    wb = waitbar(0,'Processing files (all sheets except Sheet 1)…');
    for i = 1:numel(files)
        waitbar(i/numel(files),wb,files(i).name)
        try
            fullName   = fullfile(files(i).folder,files(i).name);
            [~,sheets] = xlsfinfo(fullName);

            % Skip first sheet (index 1)
            for s = 2:numel(sheets)
                try
                    tbl = loadParticleData(fullName,'Sheet',sheets{s});

                    mdl = fitBestModel(tbl);

                    base   = erase(files(i).name,{'.xlsx','.xls'});
                    outPng = fullfile(figDir, sprintf('%s_%s.png',base,sheets{s}));
                    plotFit(tbl,mdl,outPng);

                    results = [results;
                               {sprintf('%s:%s',files(i).name,sheets{s}), ...
                                mdl.type, mdl.BIC}];
                catch innerME
                    warning("  ↳ Sheet %s skipped in %s: %s", ...
                            sheets{s}, files(i).name, innerME.message);
                end
            end
        catch ME
            warning("Skipped %s : %s",files(i).name,ME.message);
        end
    end
    close(wb)
end

% Write summary CSV (only rows gathered by the active section)

if ~isempty(results)
    resultsTbl = cell2table(results, ...
                  'VariableNames', {'File','Model','BIC'});
    writetable(resultsTbl, fullfile(rootDir,'Fit_Summary.csv'));
    disp("✔  Batch complete.  Results written to Fit_Summary.csv");
else
    disp("No results generated.  Did you disable both sections?");
end