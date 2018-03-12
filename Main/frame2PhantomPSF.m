function frame2PhantomPSF( fileNameRoot, pwv, maxFps, bIncreasingFps )

N = 2000;

if nargin < 3, maxFps = 148; end
if nargin < 4, bIncreasingFps = 0; end

numXpix = 1024; %
numYpix =  512; %

%  Load the computer phantom data
phantom_positions = tissue_pht( N, 1 );
phantom_positions(:,2) = 0; % For 2D phantom

%%  Phantom dimensions
x_size = 40/1000;            % Width of phantom [m]
% y_size = y_size/1000;        % Transverse width of phantom [m]
z_size = 20/1000;            % Height of phantom [m]
z_plus = z_size/2 + 50/1000; % Start of phantom surface [m];

%%
R = 4 / 2 / 1000; %  Radius of blood vessel [m]
Wt = R/10; % Thickness of vessel's wall
y = phantom_positions(:,2);
z = phantom_positions(:,3);
% The inside includes the wall
inside = y.^2 + (z-z_plus).^2 <= (R + Wt)^2;
outside = ~inside;
phantom_positions_out = phantom_positions(outside,:);

%% Vessel's wall
df = 30; % density factor
phantom_positions2 = tissue_pht( N*df, 1 );
y = phantom_positions2(:,2);
z = phantom_positions2(:,3);
wall = logical( ( y.^2 + (z-z_plus).^2 >=        R^2 ) .* ...
    ( y.^2 + (z-z_plus).^2 <= (R + Wt)^2 )    );

phantom_positions_wall = phantom_positions2(wall,:);

P = [ phantom_positions_out; phantom_positions_wall ];
P(:,2,:) = [];
amp = abs( randn(size(P,1),1)*4);

Q(:,1) = ( P(:,1) + 0.02 )/0.04*( numXpix - 1 ) + 1;
Q(:,2) = ( P(:,2) - 0.05 )/0.02*( numYpix - 1 ) + 1;
P = Q;

for i = 1 : size(P,1); H(i,:) = [ P(i,2)-amp(i)   P(i,2,1)+amp(i)   ]; end
for i = 1 : size(P,1); W(i,:) = [ P(i,1)-amp(i)*2 P(i,1,1)+amp(i)*2 ]; end
H = floor( H );
W = floor( W );
H( H>numYpix ) = numYpix;
H( H<1) = 1;
W( W>numXpix ) = numXpix;
W( W<1) = 1;
F = zeros( numYpix, numXpix );
for i = 1 : size(P,1)
    vH = H(i,1) : H(i,2);
    vW = W(i,1) : W(i,2);
    F( vH, vW ) = 1;
end

%%
F = F - min( F(:) );
F = F / max( F(:) );

%%
video = [];
cmppxx = 40/numXpix/10;

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
m = (p2-p1+1)/2+p1-1;

minFps = maxFps;
if bIncreasingFps, minFps = 1; end
for fps = minFps : maxFps
    disp( fps )
    nFrames = round(60/maxFps*fps);
    a = 0 : 2*pi/(nFrames+1-1) : 2*pi;
    a(end) = [];
    V = zeros( 512, 1024, nFrames );
    for t = 1 : length(a)
        parfor iC = 1 : 1024
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
            V(:,iC,t) = interp1( 1:512, Fc, zq )';
        end
    end
    
    %%
    PSF = fspecial( 'gaussian', [12 60], 20 );
    video = imfilter( V, PSF );
    
    video = video*255;
    video = round(video);
    video = uint8(video);
    
    mainDir = cd;
    mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
    dataDir = fullfile( mainDir, 'Data' );
    
    fileName = [ fileNameRoot '_fps_' num2str(fps) ];
    save( fullfile( dataDir, fileName ), 'video', 'fps', 'cmppxx', 'P', 'amp' )
end