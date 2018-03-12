clear

load NLdata

SFmm = 0.77/32; % Scaling factor sample (at 32 MHz) to mm (velcity of US in tissue 1540 m/s)
NpointsDepth = 4500;
NlinesLat  = 40;
Depth = 0 : NpointsDepth/(NpointsDepth - 1) : NpointsDepth;
AnglesRadA = -15 : 30/(NlinesLat-1) : 15;
parfor n = 1 : NlinesLat
    for m = 1 : NpointsDepth
        Xa(m,n) =  Depth(m) * sind(AnglesRadA(n));
        Ya(m,n) = -Depth(m) * cosd(AnglesRadA(n));
    end
end

AoCntr = 1574; % Estimation of the center axis of the aorta based on the markers first frame
Strt =  AoCntr - 500;
Stp  =  AoCntr + 500;

X = X/1000;
nFrames = size(X,3);
XH = zeros( length(Strt:Stp), NlinesLat, nFrames );
parfor t = 1 : nFrames
    XH(:,:,t) = abs(hilbert(X(:,:,t)));
end
figure(1)
figure('position', [ 1 41 1366 652])
clf
colordef(1, 'none')
for t = 1:nFrames
    figure(1)
    pcolor( Xa(Strt:Stp,:), Ya(Strt:Stp,:), XXH(:,:,t) );
    title(['Longitudinal plane ' num2str(t)], 'fontsize', 20 )
    nrColorbits    = 255; % Set as required
    brightenFactor = 0.5; % Set as required
    brighten( colormap( hot(nrColorbits) ), brightenFactor );
    set(gcf,'DoubleBuffer','on');
    set(gcf,'Renderer','zbuffer')
    shading('interp')
    grid on
    axis equal
    hold on
end