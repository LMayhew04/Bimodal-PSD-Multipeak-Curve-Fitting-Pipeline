function tbl = loadParticleData(fileName,varargin)
%LOADPARTICLEDATA  Import a worksheet and return a clean, normalised table
%                  with two columns:
%                     ParticleDiameter    [µm]
%                     SizeDistribution1   (area = 1 PDF)
%
%   tbl = loadParticleData(fileName)
%   tbl = loadParticleData(fileName,'Sheet','Measurement 5', ...
%                              'DiameterIdx',6,'DistIdx',8)

% ---------- User-settable defaults --------------------------------------
p = inputParser;
addParameter(p,'Sheet',2);
addParameter(p,'DiameterIdx',6,@(n)isnumeric(n)&&n>0);
addParameter(p,'DistIdx',8,@(n)isnumeric(n)&&n>0);
parse(p,varargin{:});
S = p.Results;

% ---------- Build options & select two numeric columns ------------------
opts = detectImportOptions(fileName,'Sheet',S.Sheet,...
                           'VariableNamingRule','preserve');
want = {'ParticleDiameter','SizeDistribution1'};
have = ismember(want, opts.VariableNames);

if all(have)
    opts.SelectedVariableNames = want;
else
    keepIdx = unique([S.DiameterIdx S.DistIdx]);
    if max(keepIdx) > numel(opts.VariableNames)
        error("File %s lacks columns %d & %d",fileName,S.DiameterIdx,S.DistIdx);
    end
    opts.SelectedVariableNames = opts.VariableNames(keepIdx);
end
opts = setvartype(opts,opts.SelectedVariableNames,'double');
tbl  = readtable(fileName,opts,"UseExcel",false);

if ~all(have)
    tbl.Properties.VariableNames = want;
end

% ---------- Basic clean-up ----------------------------------------------
tbl = rmmissing(tbl);
tbl = tbl(tbl.ParticleDiameter>0,:);
tbl = sortrows(tbl,'ParticleDiameter');

% ---------- NEW: normalise counts so ∫PDF dD = 1 ------------------------
area = trapz(tbl.ParticleDiameter, tbl.SizeDistribution1);
if area==0, error("%s – selected columns sum to zero.",fileName); end
tbl.SizeDistribution1 = tbl.SizeDistribution1 ./ area;
end