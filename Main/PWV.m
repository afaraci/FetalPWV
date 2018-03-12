function [ pwvX, pwvP, pwvF, d, maxPwv ] = PWV( fileName, nD, opts )

%%
if ~isfield(opts,'bDebug'),                       bDebug= 0; else,            bDebug=opts.bDebug;            end
if ~isfield(opts,'frameStep'),                 frameStep=10; else,         frameStep=opts.frameStep;         end
if ~isfield(opts,'bDebugWaves'),             bDebugWaves= 0; else,       bDebugWaves=opts.bDebugWaves;       end
if ~isfield(opts,'bPixelSelection'),     bPixelSelection= 0; else,   bPixelSelection=opts.bPixelSelection;   end
if ~isfield(opts,'bGIF'),                           bGIF= 0; else,              bGIF=opts.bGIF;              end
if ~isfield(opts,'bStabilise'),               bStabilise= 0; else,        bStabilise=opts.bStabilise;        end
if ~isfield(opts,'bRot'),                           bRot= 0; else,              bRot=opts.bRot;              end
if ~isfield(opts,'bAskRot'),                     bAskRot= 0; else,           bAskRot=opts.bAskRot;           end
if ~isfield(opts,'bIncreasingLength'), bIncreasingLength= 0; else, bIncreasingLength=opts.bIncreasingLength; end
if ~isfield(opts,'bCorrWaves'),               bCorrWaves= 0; else,        bCorrWaves=opts.bCorrWaves;        end
if ~isfield(opts,'bPlotBoxes'),               bPlotBoxes=[]; else,        bPlotBoxes=opts.bPlotBoxes;        end
if ~isfield(opts,'sequL'),                         sequL=[]; else,             sequL=opts.sequL;             end

%%
bCrop = 0;
mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
if ~( contains(fileName,'simsig') || contains(fileName,'NLdata') ...
        || contains(fileName,'phantom') )
    dataDir = fullfile( mainDir, 'Data', 'Clinical' );
else
    dataDir = fullfile( mainDir, 'Data', 'Phantom' );
end
% if contains( fileName, 'phantom' )
%     load( fullfile( dataDir, fileName ), 'dis' );
% end
fileNameRot = [ fileName '_Rot.mat' ];
cond = ~exist( 'dis', 'var' );
if cond
    if ~( contains(fileName,'simsig') || contains(fileName,'NLdata') ...
            || contains(fileName,'phantom') )
        %% Load video
        % [ video, info ] = ReadData3D;
        
        if ~exist( fullfile( dataDir, fileNameRot ), 'file' )
            video = dicomread( fileName );
            video = squeeze( video(:,:,1,:) );
            % video = mat2gray( video );
            info = dicominfo( fileName );
            cmppxx = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX;
            nFrames = info.NumberOfFrames;
            % fps = info.RecommendedDisplayFrameRate;
            fps = info.CineRate;
            save( fullfile( dataDir, fileName ), 'fps', 'cmppxx', 'nFrames' );
        else
            load( fullfile( dataDir, fileNameRot ), 'video', 'fps', 'cmppxx', 'nFrames' )
        end
    else
        if contains(fileName,'NLdata')
            load( fullfile( dataDir, 'NLdata.mat' ) )
            video = W;
            clear W V X
            nFrames = size(video,3);
        elseif contains( fileName, 'simsig' )
            load( fullfile( dataDir, fileName ), 'video', 'fps', 'cmppxx' )
            video(:,:,sequL+1:end) = [];
            nFrames = size(video,3);
        elseif contains( fileName, 'phantom' )
            if exist( fullfile( dataDir, fileNameRot ), 'file' )
                load( fullfile( dataDir, fileNameRot ), 'video', 'fps', 'cmppxx' )
                nFrames = size(video,3);
            else
                load( fullfile( dataDir, fileName ), 'video', 'fps', 'cmppxx' )
                video2 = video;
                nFrames = size(video2,3);
                for i = 0 : 2
                    video2( :, :, nFrames*i+1 : nFrames*(i+1) ) = video;
                end
                video = video2;
                video(:,:,sequL+1:end) = [];
                nFrames = size(video,3);
            end
        end
    end
    
    %% Crop video
    if bCrop
        [~,rect] = imcrop(video(:,:,1));
        video = video( rect(2):rect(2)+rect(4), ...
            rect(1):rect(1)+rect(3), : );
    end
    
    %%
    if bRot && ~exist( fullfile( dataDir, fileNameRot ), 'file' )
        imshow(video(:,:,1));
        if bAskRot
            str = input('Would you like to rotate the image?(Y/N): ','s');
        else
            str = 'N';
        end
        if str == 'Y'
            [ c, r, ~ ] = impixel( video(:,:,1) );
            tang = ( r(1) - r(2) ) / ( c(1) - c(2) );
            angle = atand(tang);
            warning('off')
            fprintf('\n');
            fprintf('Progress:                ');
            parfor iFrame = 1 : nFrames
                fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' );
                fprintf( '%3.0f%% of frames.', floor(iFrame/nFrames*100) );
                Rvideo(:,:,iFrame) = imrotate(video(:,:,iFrame),angle); %#ok<AGROW>
                % Rvideo = imrotate3(video,angle,[ 0 0 1]); % It might take
                % too long!
            end
            video = Rvideo;
        end
        save( fullfile( dataDir, fileNameRot ), 'video', 'cmppxx', 'nFrames', 'fps' )
    end
    close all
    
    if bStabilise
        video = stabiliseVideo(video,cmppxx);
    end
    
    warning('off')
    fprintf('\n');
    fprintf('Progress:                ');
    for iFrame = 1 : nFrames
        fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' );
        fprintf( '%3.0f%% of frames.', floor(iFrame/nFrames*100) );
        vFrame = imadjust( wiener2(video(:,:,iFrame),[5 5]) );
        if iFrame == 1
            % Pixel selection
            if bPixelSelection
                [ c, r, ~ ] = impixel( vFrame );
                save( fullfile( dataDir, fileName ),'c','r')
            elseif exist( fullfile( dataDir, [ fileName '.mat' ] ),'file')
                load( fullfile( dataDir, fileName ),'c','r')
            end
            c1 = min(c);
            cE = max(c);
            
            maxPwv0 = ( cE - c1 - 1 ) * cmppxx * 10 * fps / 1000; % [m/s]
            fprintf('\n')
            fprintf('Maximum initial PWV = %2.2f m/s \n\n',maxPwv0 );
            
            save( fullfile( dataDir, fileName ), 'maxPwv0', '-append' );
            pause(2);
            fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' );
            fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' );
            d = unique(round( c1 : (cE-c1)/(nD-1) : cE ))';
            nD = length( d );
            c = zeros(nD*2,1);
            c(1:2:end) = d;
            c(2:2:end) = d;
            nP = nD*2; % Number of points
            dis = zeros(nD,nFrames);
            cM  = zeros(nD,nFrames);
            rM  = zeros(nD,nFrames);
            mmwH = 4; % Window height [mm]
            mmwW = 0.125; % Window width [mm]
            AborH = round(mmwH/(cmppxx*10)/2);
            AborW = round(mmwW/(cmppxx*10)/2);
            A = zeros( AborH*2+1, AborW*2+1, nP );
            r1 = r(1);
            rE = r(end);
            r(1:2:nD*2) = r1;
            r(2:2:nD*2) = rE;
            for iP = 1 : nP
                A(:,:,iP) = vFrame( r(iP)-AborH:r(iP)+AborH, c(iP)-AborW:c(iP)+AborW );
            end
            % Move each point to the edge of the vessel's walls
            r = move2edge( r, AborH, nD, A );
            r = reshape(r,2,length(r)/2);
            TF1 = is_outlier( r(1,:), 2, 'median' );
            r(1,:) = interp1( d(~TF1)', r(1,~TF1), d', 'nearest', 'extrap' );
            TF2 = is_outlier( r(2,:), 2, 'median' );
            r(2,:) = interp1( d(~TF2)', r(2,~TF2), d', 'nearest', 'extrap' );
            r = reshape(r,length(r)*2,1);
        else
            for iP = 1 : nP
                A(:,:,iP) = vFrame( r(iP)-AborH:r(iP)+AborH, c(iP)-AborW:c(iP)+AborW );
            end
            % Move each point to the edge of the vessel's walls
            r = move2edge( r, AborH, nD, A );
            r = reshape(r,2,length(r)/2);
            TF1 = is_outlier( r(1,:), 2, 'median' );
            r(1,TF1) = nan;
            tmp = ...
                ( fillmissing([ rU(iFrame-1, : ); r(1, : )]', 'nearest')' + ...
                fillmissing([ rU(iFrame-1, : ); r(1, : )],  'nearest')  )/2;
            r(1,:) = tmp(2,:);
            spanmm = 0.5; % [mm]
            span = round( ( spanmm/ ( mean(diff(c(1:2:end)))*cmppxx*10 ) -1 )/ 2 ) * 2 + 1;
            tmp = imfilter([ rU(iFrame-1, : ); r(1, : )],fspecial('average', [2 span]),'replicate');
            r(1,:) = tmp(1,:);
            TF2 = is_outlier( r(2,:), 2, 'median' );
            r(2,TF2) = nan;
            tmp = ...
                ( fillmissing([ rL(iFrame-1, : ); r(2, : )]', 'nearest')' + ...
                fillmissing([ rL(iFrame-1, : ); r(2, : )],  'nearest')  )/2;
            r(2,:) = tmp(2,:);
            tmp = imfilter([ rL(iFrame-1, : ); r(2, : )],fspecial('average', [2 span]),'replicate');
            r(2,:) = tmp(1,:);
            r = reshape(r,length(r)*2,1);
        end
        %[r, c] = smoothEdge(r,c,nD);
        P = [c r];
        
        rU(iFrame,:) = r(1:2:end);
        rL(iFrame,:) = r(2:2:end);
        
        % Distance between points
        dis(:,iFrame) = abs( r(1:2:end) - r(2:2:end) );
        cM(:,iFrame) = ( P(1:2:end,1) + P(2:2:end,1) ) / 2;
        rM(:,iFrame) = ( P(1:2:end,2) + P(2:2:end,2) ) / 2;
        if bDebug && iFrame == round( (iFrame-1)/frameStep )*frameStep + 1
            imshow( vFrame, 'InitialMagnification', 100 );
            hold on
            plot(c,r,'ro','MarkerFaceColor','r')
            plot(cM(:,iFrame),rM(:,iFrame),'g-','Marker','x')
            for iD = 1 : nD
                i = iD*2-1;
                plot(P(i:i+1,1),P(i:i+1,2),'b-');
            end
            if bPlotBoxes
                for iP = 1 : nP
                    crnrs = [ c(iP)-AborW r(iP)-AborH
                        c(iP)+AborW r(iP)-AborH
                        c(iP)+AborW r(iP)+AborH
                        c(iP)-AborW r(iP)+AborH ];
                    crnrs(end+1,:) = crnrs(1,:);
                    line(crnrs(:,1),crnrs(:,2),'LineWidth',2,'Color','r')
                end
            end
            drawnow
            if bGIF
                % Create the gif file
                gifname = 'animation.gif';
                frame = getframe(1);
                im = frame2im(frame);
                [im,map] = rgb2ind(im,256);
                if iFrame == 1
                    imwrite(im,map,gifname,'gif','LoopCount',Inf,'DelayTime', 0.05);
                else
                    imwrite(im,map,gifname,'gif','WriteMode','append','DelayTime', 0.05);
                end
            end
        end
    end
    save( fullfile( dataDir, fileName ), 'fps', 'cmppxx', 'dis', 'd', 'c1', 'cE', 'maxPwv0', '-append' );
else
    fps = 0; cmppxx = 0; c1 = 0; cE = 0;
    load( fullfile( dataDir, fileName ), 'd', 'dis', 'fps', 'cmppxx','c1', 'cE' );
    maxPwv0 = ( cE - c1 - 1 ) * cmppxx * 10 * fps / 1000; % [m/s]
    disp( maxPwv0 );
    nFrames = size( dis, 2 );
end

%%
t = ( 1 : nFrames ) / fps;
dis = dis';

%% Plot original waveforms
if bDebugWaves
    figure
    subplot(2,1,1)
    plot( t, dis*cmppxx*10 )
    title('Diameter vs. Time');
    ylabel('Diameter (mm)');
    xlabel('Time (s)');
    axis tight
    drawnow
end

%% Filter waves
nW = size( dis, 2 );
warning('off')
fprintf('\n');
fprintf('Progress:                ');
ord = floor( size( dis, 1 ) / 4 );
StBfF = 2;
for iW = 1 : size( dis, 2 )
    fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' );
    fprintf( '%3.0f%% of waves. ', floor(iW/nW*100) );
    dis(:,iW) = filterWave( dis(:,iW), StBfF, ord );
    %% Plot filtered waveforms
    if bDebugWaves
        subplot(2,1,2)
        plot( t, dis(:,iW)*cmppxx*10 )
        title('Diameter vs. Time - Filtered waves');
        ylabel('Diameter (mm)');
        xlabel('Time (s)');
        axis tight
        hold on
        drawnow
    end
end

%% Detrending of waveforms
detrendType = 'constant';
dis = detrend(dis,detrendType);

%% Standardisation of waveforms
dis = zscore(dis);

%% Selection of correlated waveforms
if bCorrWaves
    Zdis = zscore(dis);
    Rtst1 = corr(Zdis);
    Rset = 0.67; % 0.985 for testing on simsig and phantom
    tmp = Rtst1 > Rset;
    Rtst2 = tmp.* Rtst1;
    Rtst3 = sum(Rtst2);
    [~, Index] = max(Rtst3);
    Rtst = tmp(Index,:); % select the lines exhibiting the most and highest correlations
    LineSelect = find( Rtst == 1 );
    if length( LineSelect ) > 1
        dis = dis(:,LineSelect);
        d = d(LineSelect);
    end
end

%% Save final waves
disF = dis;
dF = d;
save( fullfile( dataDir, fileName ), 'dF', 'disF', '-append' );

%% Plot of final waveforms
if bDebugWaves
    figure
    plot( t, dis*cmppxx*10 )
    title('Diameter vs. Time');
    ylabel('Diameter (mm)');
    xlabel('Time (s)');
    axis tight
    drawnow
end

%% Cross correlation
bDebug = 1;
[ pwvX, maxPwv ] = crossCorr( fileName, d, dis, cmppxx, fps, bIncreasingLength, bDebug );
[ pwvP, pwvF ] = ttp_ttf( dis, fps, cmppxx, d );

%% Functions
    function r = move2edge( r, AborH, nD, A )
        v = ( 1 : AborH*2 + 1 )';
        bB = round(AborH/4);
        bT = round(AborH/4);
        w = ( v >= AborH + 1 - bB ) .* ( v <= AborH + 1 + bT );
        smf = round(round(AborH/4)/2)*2 - 1;
        for j = 0 : 1
            for k = 1 : nD
                kk = k*2-1 + j;
                tmpA = A(:,:,kk);
                tmpSig = smooth( mean(tmpA,2), smf );
                m = max( diff(tmpSig,1).*w(2:end)*(2*j-1) );
                idx = floor( mean( find( diff(tmpSig,1).*w(2:end)*(2*j-1) == m ) ) );
                idx = idx + 1 - (2*j-1) ;
                r(kk) = r(kk) + idx - AborH - 1;
            end
        end
    end

    function [r, c] = smoothEdge(r,c,nD)
        %         r(1:2:end) = round(smooth(r(1:2:end),3,'moving'));
        %         r(2:2:end) = round(smooth(r(2:2:end),3,'moving'));
        p1 = polyfit(c(1:2:end),r(1:2:end),5);
        p2 = polyfit(c(2:2:end),r(2:2:end),5);
        r(1:2:end) = polyval( p1, c(1:2:end) );
        r(2:2:end) = polyval( p2, c(2:2:end) );
        for k = 1 : nD
            j = k*2-1;
            c(j) = mean(c(j:j+1));
            c(j+1) = c(j);
        end
    end

    function TF = is_outlier( A, f, method )
        if nargin < 2
            f = 3;
        end
        if nargin < 3
            method = 'median';
        end
        K = 1.4826;
        switch method
            case 'median'
                AD = abs( A - median( A ) ); % Absolute deviation
                MAD = median( AD ); % Median of absolute deviation
            case 'mean'
                AD = abs( A - mean( A ) ); % Absolute deviation
                MAD = mean( AD ); % Mean of absolute deviation
        end
        TF = AD > f*K*MAD;
    end
end