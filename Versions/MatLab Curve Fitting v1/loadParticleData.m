function tbl = loadParticleData(fileName,varargin)
% Load Particle Data: Import any workbook and return a table that
% with the two columns (even if nonexistent in raw spreadsheet):
%   ParticleDiameter
%   SizeDistribution1
%
%   tbl = loadParticleData(fileName)
%   tbl = loadParticleData(fileName,'Sheet','Measurement 5', ...
%                              'DiameterIdx',6,'DistIdx',8)
%
% Optional name-value pairs
%   'Sheet'        – target worksheet (default = 2nd sheet)
%   'DiameterIdx'  – column index to take as particle diameter  (default = 6)
%   'DistIdx'      – column index to take as size distribution   (default = 8)

%% Defaults and user overrides
p = inputParser;
addParameter(p,'Sheet',2);
addParameter(p,'DiameterIdx',6,@(n) isnumeric(n) && n>0);
addParameter(p,'DistIdx',8,@(n) isnumeric(n) && n>0);
parse(p,varargin{:});
S = p.Results;

%% Set up import

% 1. Build import options (auto-detects types, ranges, etc.)
opts = detectImportOptions(fileName,'Sheet',S.Sheet,...
                           'VariableNamingRule','preserve');

% 2. Make sure we grab exactly two numeric columns.
%    If the workbook *already* has the right names, keep them.
haveNames = ismember({'ParticleDiameter','SizeDistribution1'},opts.VariableNames);

if all(haveNames)
    opts.SelectedVariableNames = {'ParticleDiameter','SizeDistribution1'};

else
    % -- Fallback: pull by column index, then rename like the wizard did
    keepIdx = unique([S.DiameterIdx, S.DistIdx]);
    if max(keepIdx) > numel(opts.VariableNames)
        error("Workbook %s does not have columns %d and %d.",...
              fileName,S.DiameterIdx,S.DistIdx);
    end
    opts.SelectedVariableNames = opts.VariableNames(keepIdx);
end

opts = setvartype(opts,opts.SelectedVariableNames,'double');   % force numeric
tbl  = readtable(fileName,opts,"UseExcel",false);

% 3. Rename if necessary so downstream code is 100 % stable
if ~all(haveNames)
    tbl.Properties.VariableNames = {'ParticleDiameter','SizeDistribution1'};
end

%% 4. Basic clean-up
% (sort, force positive values, etc.)
tbl  = rmmissing(tbl);
tbl  = tbl(tbl.ParticleDiameter>0,:);
tbl  = sortrows(tbl,'ParticleDiameter');
end
