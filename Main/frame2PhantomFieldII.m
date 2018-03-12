function frame2PhantomFieldII( fileName, pwv )

mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
dataDir = fullfile( mainDir, 'Data' );

video = [];
load( fullfile( dataDir, fileName ),'video','fps','cmppxx')
F = double( video(:,:,1) );

hi = size( F, 1 );
him = hi/2;
a = him-hi*25/2^6;
hia = him-a+1;
hib = him+a-1;
[ ~, Mloc1 ] = max( diff(mean( F(hia:him,:),2), 2 ) );
p1 = Mloc1 + 2 + hia - 1;
[ ~, Mloc2 ] = max( diff(mean( F(him+1:hib+1,:),2), 2 ) );
p2 = Mloc2 + 2 + him;
if mod( (p2-1)-p1, 2), p2 = p2-1; end

nFrames = 60;
m = (p2-p1+1)/2+p1-1;
a = 0 : 2*pi/(nFrames+1-1) : 2*pi;
a(end) = [];

vq = zeros( 512, 1024, nFrames );
for t = 1 : length(a)
    t
    for iC = 1 : 1024
        b = (iC-1)*cmppxx*10/(pwv*1000)*fps;
        b = (b-1)/nFrames*(2*pi);
        np1 = round( p1-sin(a-b)*5 );
        np2 = round( p2+sin(a-b)*5 );
        x1q =  1 : (p1-1)/(np1(t)-1) : p1;
        x2q = p1 : (m-p1)/((m-np1(t)+1)-1) : m;
        xq = [x1q x2q(2:end)];
        y1q =  p2 : (512-p2)/(512-np2(t)) : 512;
        y2q = m+1 : (p2-m-1)/((512-m-(512-np2(t)))-1) : p2;
        yq = [ y2q(1:end-1) y1q ];
        zq = [ xq yq ]';
        Fc = F(:,iC);
        vq(:,iC,t) = interp1( 1:512, Fc, zq )';
    end
end
video = uint8(vq);
fileName = [ fileName 'C' ];
save( fullfile( dataDir, fileName ), 'video', 'fps', 'cmppxx' )