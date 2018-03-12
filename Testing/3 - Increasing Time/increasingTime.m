function [ pwv, sequL ] = increasingTime( fileName, minSequL, maxSequL )

if nargin < 2, minSequL = 51; end
if nargin < 3, maxSequL = 363; end

nPoints = 338;
load(fileName,'video','fps','cmppxx')
maxC = calculateC( cmppxx, size(video,2), nPoints );
opts.inputC = maxC;
opts.bPixelSelection = 0;
opts.bXCorr = 1;
opts.bTimeSmooth = 0;
opts.bDebug = 0; opts.frameStep = 10;
opts.bDebugWaves = 0;
opts.bDebugResults = 0;
opts.aortaLength = [];
opts.bPlotBoxes = 0;
sequL = minSequL : maxSequL;
nTimes = length( sequL );
pwv = zeros(1,nTimes);
for iPWV = 1 : nTimes
    disp(iPWV)
    opts.sequL = sequL( iPWV );
    pwv(iPWV) = PWV( fileName, nPoints, 1, opts );
end