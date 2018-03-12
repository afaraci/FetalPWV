function pwv = increasingPoints( fileName )

mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
dataDir = fullfile( mainDir, 'Data' );
load( fullfile( dataDir, fileName ),'dis');
maxNpoints = size(dis,1);

pwv = zeros( 1, maxNpoints );
opts = [];
minNpoints = 2;
parfor nPoints = minNpoints : maxNpoints
    display(nPoints)
    pwv(nPoints) =  PWV( fileName, nPoints, opts );   
end

pwv(1) = [];