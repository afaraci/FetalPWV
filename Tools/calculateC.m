function c = calculateC( fileName, nD, bPixelSelection )

if nargin < 3
    bPixelSelection = 0;
end

mainDir = cd;
dataDir = fullfile( mainDir, 'Data' );

video = dicomread( fileName );
video = squeeze( video(:,:,1,:) );
video = mat2gray( video );
info = dicominfo( fileName );
cmppxx = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX;
vFrame = video(:,:,1);

if bPixelSelection
    c = impixel( vFrame );
    mkdir(fullfile(dataDir,'Input','TrackPoints'));
    save(fullfile(dataDir,'Input','TrackPoints',[ fileName '.mat'] ),'c','r')
elseif exist(fullfile(dataDir,'Input','TrackPoints',[ fileName '.mat'] ),'file')
    load(fullfile(dataDir,'Input','TrackPoints',[ fileName '.mat'] ),'c','r')
end
c1 = c(1);
cE = c(end);
for iD = 1 : nD
    i = iD*2-1;
    c(i) = round( (cE-c1)/(nD-1)*(iD-1)+c1 );
    c(i+1) = c(i);
end