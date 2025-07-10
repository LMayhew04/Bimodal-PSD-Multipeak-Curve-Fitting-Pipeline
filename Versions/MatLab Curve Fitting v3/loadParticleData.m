function tbl = loadParticleData(fileName, varargin)
%LOADPARTICLEDATA  Import a worksheet and return a clean, normalized table.
%
%   tbl = loadParticleData(fileName)
%   tbl = loadParticleData(fileName, 'Sheet', 2, 'DiameterIdx', 6, 'DistIdx', 8)
%   tbl = loadParticleData(..., 'CustomNames', {'MyDiameter', 'MyDistribution'})
%
%   This function returns a table with two columns:
%     - ParticleDiameter (µm)
%     - SizeDistribution1 (PDF normalized to area = 1)
%
%   Optional name-value pairs:
%     'Sheet'        – worksheet number or name (default: 2)
%     'DiameterIdx'  – index of diameter column (default: 6)
%     'DistIdx'      – index of distribution column (default: 8)
%     'CustomNames'  – cell array with 2 strings: desired column names

% Parse input options
p = inputParser;
addParameter(p, 'Sheet', 2);
addParameter(p, 'DiameterIdx', 6, @(n) isnumeric(n) && n > 0);
addParameter(p, 'DistIdx', 8, @(n) isnumeric(n) && n > 0);
addParameter(p, 'CustomNames', {'ParticleDiameter', 'SizeDistribution1'});
parse(p, varargin{:});
S = p.Results;

% Try to find named columns, fallback to index
opts = detectImportOptions(fileName, 'Sheet', S.Sheet, ...
    'VariableNamingRule', 'preserve');

haveNames = ismember(S.CustomNames, opts.VariableNames);

if all(haveNames)
    opts.SelectedVariableNames = S.CustomNames;
    raw = readtable(fileName, opts);
    raw.Properties.VariableNames = {'ParticleDiameter', 'SizeDistribution1'};
else
    idx = unique([S.DiameterIdx S.DistIdx]);
    if max(idx) > numel(opts.VariableNames)
        error("File %s lacks columns %d & %d", fileName, S.DiameterIdx, S.DistIdx);
    end
    opts.SelectedVariableNames = opts.VariableNames(idx);
    raw = readtable(fileName, opts);
    raw.Properties.VariableNames = {'ParticleDiameter', 'SizeDistribution1'};
end

% Clean data
tbl = raw(:, {'ParticleDiameter', 'SizeDistribution1'});
tbl = rmmissing(tbl);
tbl = sortrows(tbl, 'ParticleDiameter');

% Normalize area under PDF to 1
x = tbl.ParticleDiameter;
y = tbl.SizeDistribution1;
area = trapz(x, y);
tbl.SizeDistribution1 = y / area;
end
