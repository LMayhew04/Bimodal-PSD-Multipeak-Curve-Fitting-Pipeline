function tbl = loadParticleData(filename, sheet)
%LOADPARTICLEDATA  Loads and preprocesses a particle size sheet.

data = readtable(filename, 'Sheet', sheet);
headers = lower(data.Properties.VariableNames);

idxDiam = find(contains(headers, 'diam'), 1);
idxDist = find(contains(headers, 'dist'), 1);

if isempty(idxDiam) || isempty(idxDist)
    % Fallback: use columns 1 and 2
    idxDiam = 1; idxDist = 2;
end

ParticleDiameter = data{:, idxDiam};
SizeDistribution1 = data{:, idxDist};

% Remove NaN/empty rows
good = ~isnan(ParticleDiameter) & ~isnan(SizeDistribution1);
ParticleDiameter = ParticleDiameter(good);
SizeDistribution1 = SizeDistribution1(good);

tbl = table(ParticleDiameter, SizeDistribution1);

end