function createPhantom( nameSufx )

bVideo = 0;
mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
dataDir = fullfile( mainDir, 'Data' );

fileName = [ 'phantom' nameSufx ];

%%
pwv = -3;    % Pulse Wave Velocity [m/s]
fps = 148;   % Frames per second [Hz]
nLines = 50;
[ rf_lines, tstart_lines ] = generatePWVphantomImage( pwv, fps, nLines );
% save( fullfile( dataDir, fileName ), 'rf_lines', 'tstart_lines');
nLines  = size( rf_lines, 1 );
nFrames = size( rf_lines, 2 );
 
%% Parameteres
fs = 100e6; % Sampling frequency [Hz]
c = 1540;   % Speed of sound [m/s]
min_sample = 0;
image_width = 40;               % Size of image sector [mm]
image_width = image_width/1000; % Size of image sector [m]
dx = image_width/nLines;        % Increment for image  [m/line]

%%
for iFrame = 1 : nFrames
    for iLine = 1 : nLines
        rf_data = cell2mat( rf_lines(iLine,iFrame) );
        tstart = tstart_lines( iLine, iFrame );
        if tstart > 0
            rf = [ zeros( floor(tstart*fs)-min_sample, 1 ); rf_data ];
        elseif tstart == 0
            rf = rf_data( abs(tstart*fs)+1 : length(rf_data) );
        end
        rfFrame( 1 : length(rf), iLine ) = rf;
    end
    while iFrame > 1 && size(rfFrame,1)> size(frames,1)
        rfFrame(1,:) = [];
    end
    frames(:,:,iFrame) = abs( hilbert( rfFrame ) );
end

%%
frames = frames - min( frames(:) );
frames = frames / max( frames(:) );
frames = log( frames + 0.01 );
[ lLines, nLines, nFrames ] = size( frames );

%%
xSpace = (   1:nLines )*dx*1000; % [mm]
ySpace = ( ( 1:lLines )/fs + min_sample/fs )*c/2*1000; % [mm]

%%
limT = 51; % [mm]
limB = 71; % [mm]
[~,idxT] = min(abs(ySpace - limT));
[~,idxB] = min(abs(ySpace - limB));

ySpace(idxB+1: end)= [];
ySpace(1: idxT-1)= [];

frames(idxB+1: end,:,:)= [];
frames(1: idxT-1,:,:)= [];

%%
for iFrame = 1 : nFrames
    framesRsz(:,:,iFrame) = imresize(frames(:,:,iFrame),[512 1024]);
end
frames = framesRsz; clear framesRsz
xSpace = imresize( xSpace, [1 1024], 'nearest' );
ySpace = imresize( ySpace, [1  512], 'nearest' );

%%
frames = frames - min( frames(:) );
frames = frames / max( frames(:) );
frames = frames*255;
frames = round(frames);
frames = uint8(frames);

%%
video = frames;
cmppxx = ( xSpace(end) - xSpace(1) ) / 1024 / 10; % [cm]
cmppxy = ( ySpace(end) - ySpace(1) ) /  512 / 10; % [cm]
save( fullfile( dataDir, fileName ),'video','fps','cmppxx')

%% Create Video
if bVideo
    xSpace = (0:1023)*cmppxx*10; % [mm]
    ySpace = (0: 511)*cmppxy*10; % [mm]
    writerObj = ... 
        VideoWriter( fullfile( mainDir, 'Resources', [ fileName '.avi'] ), ... 
                     'Uncompressed AVI' );
    writerObj.FrameRate = fps;
    open(writerObj);
    
    figure('units','pixels','position', [ 1 1 1024 512 ], 'resize', 'off', 'toolbar', 'none' )
    axis ij
    axis equal
    axis tight
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    nFrames = size( video, 3 );
    for iFrame = 1 : nFrames
        imagesc( xSpace, ySpace, video(:,:,iFrame) )
        colormap( gray(256) )
        drawnow
        frame = getframe;
        writeVideo(writerObj,frame);
    end
    
    close(writerObj);
    close
end