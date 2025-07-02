% ===============================
% User Configuration Section
% ===============================
inputFile = '*';  % '*' or '' to run all .xlsx files in 'data/', or set filename for single file
sheetList = {2, 4, 6, 8};        % Can be sheet names or indices (cell or numeric array)

diamCol = 6;         % F column: typically particle diameter
distCol = 8;         % H column: e.g. volume-weighted distribution
dataStartRow = 8;    % Data starts in row 8 (header above)

outFolder = 'Figures/';
summaryFile = 'Fit_Summary.csv';

fitOpts.numComp = 3;
fitOpts.type = 'mix';      % 'mix', 'logn', or 'norm'
fitOpts.minWidth = 0.03;   % Min allowed width (as fraction of diameter range)
fitOpts.maxTries = 20;
fitOpts.peakTol = 0.07;    % Peak match tolerance as fraction of diameter range

% ===============================
% File Discovery & Main Batch Loop
% ===============================
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% --- Determine which Excel file(s) to process ---
if isempty(inputFile) || strcmp(inputFile, '*')
    dataFiles = dir('data/*.xlsx');
    % Filter out hidden/temporary files (like ~$Book1.xlsx)
    dataFiles = dataFiles(arrayfun(@(f) ~startsWith(f.name,'~$'), dataFiles));
    if isempty(dataFiles)
        error('No .xlsx files found in the data/ folder.');
    end
    fileList = fullfile({dataFiles.folder}, {dataFiles.name});
else
    fileList = {inputFile};
end

summary = {};
summaryHeader = {'File','Sheet','RMSE','MissedPeakPenalty','WidthPenalty','NumMissedPeaks','Params'};

for f = 1:length(fileList)
    thisFile = fileList{f};
    for s = 1:length(sheetList)
        % Robust sheet reference (handles cell array or numeric)
        if iscell(sheetList)
            thisSheet = sheetList{s};
        else
            thisSheet = sheetList(s);
        end

        try
            tbl = loadParticleDataCustom(thisFile, thisSheet, diamCol, distCol, dataStartRow);
        catch err
            warning('Failed on %s:%s - %s', thisFile, toString(thisSheet), err.message);
            continue
        end

        % Use only nonzero-distribution rows for all fitting and plotting
        good = tbl.SizeDistribution1 > 0;
        if sum(good) < 3
            warning('Not enough valid rows in %s, sheet %s. Skipping.', thisFile, toString(thisSheet));
            continue
        end

        % --- v4-style preprocessing: sorted, normalized, interpolated PDF ---
        diam = tbl.ParticleDiameter(good);
        pdf_data = tbl.SizeDistribution1(good);
        [diam, idx] = sort(diam);  % Ensure monotonically increasing
        pdf_data = pdf_data(idx);
        pdf_data = pdf_data / trapz(diam, pdf_data); % Area normalize

        % Fit the model to the PDF
        model = fitBestModel(diam, pdf_data, fitOpts);

        % Plotting
        figBase = sprintf('%s_%s_fit.png', getFileBase(thisFile), toString(thisSheet));
        outPng = fullfile(outFolder, figBase);
        plotFit(diam, pdf_data, model, outPng);

        % Diagnostics (missed peaks, width penalties)
        peak_prom = max(pdf_data)*0.05;
        dataLocs = simplePeaks(diam, pdf_data, peak_prom);
        xq = linspace(min(diam), max(diam), 400);
        y_model = model.pdf(xq);
        modelLocs = simplePeaks(xq, y_model, max(y_model)*0.04);

        missed = 0;
        for d = dataLocs(:)'
            if all(abs(modelLocs - d) > fitOpts.peakTol*range(diam))
                missed = missed + 1;
            end
        end
        widthPenalty = sum([model.components.sigma] < fitOpts.minWidth);
        rmse = sqrt(mean((interp1(xq, y_model, diam, 'linear', 'extrap')-pdf_data).^2));
        summary(end+1,:) = {getFileBase(thisFile), toString(thisSheet), rmse, missed, widthPenalty, missed, mat2str(model.params,4)};
    end
end

fid = fopen(summaryFile, 'w');
fprintf(fid, '%s,', summaryHeader{1:end-1}); fprintf(fid, '%s\n', summaryHeader{end});
for i=1:size(summary,1)
    for j=1:size(summary,2)-1
        fprintf(fid, '%s,', num2str(summary{i,j}));
    end
    fprintf(fid, '%s\n', summary{i,end});
end
fclose(fid);

% ===============================
% Helper Functions
% ===============================
function tbl = loadParticleDataCustom(inputFile, sheet, diamCol, distCol, startRow)
data = readmatrix(inputFile, 'Sheet', sheet, 'Range', sprintf('A%d:Z10000', startRow));
if size(data,2) < max(diamCol,distCol)
    error('Selected columns out of bounds');
end
diam = data(:,diamCol-1+1);  % Adjust for matrix indexing (A=1)
dist = data(:,distCol-1+1);
valid = isfinite(diam) & isfinite(dist) & ~isnan(diam) & ~isnan(dist) & (dist > 0);
diam = diam(valid);
dist = dist(valid);
tbl = table(diam, dist, 'VariableNames', {'ParticleDiameter','SizeDistribution1'});
end

function s = toString(val)
if isnumeric(val)
    s = num2str(val);
else
    s = char(val);
end
end

function base = getFileBase(filename)
[~, base, ~] = fileparts(filename);
end

function peakLocs = simplePeaks(x, y, minProm)
if nargin < 3, minProm = 0; end
peakLocs = [];
for k = 2:length(y)-1
    if y(k) > y(k-1) && y(k) > y(k+1) && y(k) > minProm
        peakLocs(end+1) = x(k); %#ok<AGROW>
    end
end
end