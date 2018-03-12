function video = stabiliseVideo(video,cmppxx,bDebug,bGif)

if nargin < 3; bDebug = 0; end
if nargin < 4;   bGif = 0; end

nFrames = size(video,3);

for iFrame = 1 : nFrames
    fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' );
    fprintf( '%3.0f%% of frames.', floor(iFrame/nFrames*100) );
    vFrame = video(:,:,iFrame);
    if iFrame == 1
        % Pixel selection
        [ c, r, ~ ] = impixel( vFrame );
        nP = length(c);
        P = zeros(nP,2,nFrames);
        % Correlation matrix
        Mcor = ones(nFrames,nP);
        mmwH = 5; % Dimension of window hight in mm
        mmwW = 5; % Dimension of window width in mm
        %         AborH = 3;
        %         AborW = 3;
        AborH = round(mmwH/(cmppxx*10)/2);
        AborW = round(mmwW/(cmppxx*10)/2);
        BborH = AborH + 4;
        BborW = AborW + 4;
        A = zeros(AborH*2+1,AborW*2+1,nP); % square box
        for iP = 1 : nP
            A(:,:,iP) = ...
                vFrame( r(iP)-AborH:r(iP)+AborH, c(iP)-AborW:c(iP)+AborW );
        end
    else
        for iP = 1 : nP
            B = vFrame( r(iP)-BborH:r(iP)+BborH, c(iP)-BborW:c(iP)+BborW );
            C = xcorr2n(A(:,:,iP),B);
            [ Mcor(iFrame,iP), I ] = max( C(:) );
            [ I_row, I_col ] = ind2sub( size(C), I );
            r(iP) = r(iP) + I_row - ( BborH-AborH + 1 );
            c(iP) = c(iP) + I_col - ( BborW-AborW + 1 );
        end
        for iP = 1 : nP
            A(:,:,iP) = ...
                vFrame(r(iP)-AborH:r(iP)+AborH,c(iP)-AborW:c(iP)+AborW);
        end
    end
    
    cc(:,iFrame) = c;
    rr(:,iFrame) = r;
    
    if bDebug || iFrame == round( (iFrame-1) /20)*20 + 1
        imshow( vFrame, 'InitialMagnification', 100 );
        hold on
        plot(c,r,'ro','MarkerFaceColor','r')
        drawnow
        if bGif
            % Create the gif file
            gifname = 'animation.gif';
            frame = getframe(1);
            im = frame2im(frame);
            [im,map] = rgb2ind(im,256);
            if iFrame == 1;
                imwrite(im,map,gifname,'gif','LoopCount',Inf,'DelayTime', 0.05);
            else
                imwrite(im,map,gifname,'gif','WriteMode','append','DelayTime', 0.05);
            end
        end
    end
end

ccSelect = selectCorrWaves(cc(2:end,:),0.9);
rrSelect = selectCorrWaves(rr(2:end,:),0.9);
x = zeros(1,nP-1);
y = zeros(1,nP-1);
x(ccSelect) = 1;
y(rrSelect) = 1;
z = x+y;
select = find(z>0);
select = select + 1;
cc = [ cc(1,:); cc(select,:)];
rr = [ rr(1,:); rr(select,:)];
c = round( mean(cc,1) );
r = round( mean(rr,1) );
cI =  (max(c) - c(1)); % maximum increase
cD = -(min(c) - c(1)); % maximum decrease
rI =  (max(r) - r(1)); % maximum increase
rD = -(min(r) - r(1)); % maximum decrease
cCh = c-c(1);
rCh = r-r(1);
videoTmp = zeros(size(video,1)+rI+rD,size(video,2)+cI+cD,size(video,3));
videoTmp(rI+1:rI+size(video,1),cI+1:cI+size(video,2),:) = video;
for iFrame = 2 : nFrames
    to = rI - rCh(iFrame) + 1;
    bo = to + size(video,1) - 1;
    le = cI - cCh(iFrame) + 1;
    ri = le + size(video,2) - 1;
    videoTmp( to:bo, le:ri, iFrame ) = video(:,:,iFrame);
end
video = videoTmp;
bDebug = 1;
if bDebug
    for iFrame = 1 : 6 : nFrames
        imshow( video( :, :, iFrame ), 'InitialMagnification', 100 );
        drawnow
    end
end

    function  waveSelect = selectCorrWaves(dd,Rset)
        %% Selection of correlated waveforms
        if nargin < 2
            Rset = 0.9;
        end
        Z = zscore(dd');
        Rtst1 = corr( Z );
        tmp = Rtst1 > Rset;
        Rtst2 = tmp.* Rtst1;
        Rtst3 = sum(Rtst2);
        [~, Index] = max(Rtst3);
        Rtst = tmp(Index,:); % select the lines exhibiting the most and highest correlations
        waveSelect = find( Rtst == 1 );
    end

end