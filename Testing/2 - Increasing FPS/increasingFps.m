function [ pwv, maxPwv ] = increasingFps( fileName, nPoints )

mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
dataDir = fullfile( mainDir, 'Data' );

video = [];
load( fullfile(dataDir,fileName), 'video', 'fps', 'cmppxx' )
video(:,:,end+1) = video(:,:,1) ;
nFrames = size( video, 3 );
nFramesNew = 1 : nFrames-1;
for i = 42 : length( nFramesNew )
    if ~exist( fullfile(dataDir,[ fileName '_' num2str(nFramesNew(i)) 'frames.mat']),'file')
        disp(i)
        if ~exist('X','var')
            [X,Y,Z] = meshgrid( 1:1024, 1:512, 1:nFrames);
        end
        [Xq,Yq,Zq] = meshgrid( 1:1024, 1:512, ...
            1 : (nFrames-1)/(nFramesNew(i)) : nFrames);
        videoOld = video;
        video = uint8( interp3(X,Y,Z,double(video),Xq,Yq,Zq) ); %#ok<NASGU>
        video(:,:,end) = [];
        fpsOld = fps;
        fps = fps/(nFrames-1)*nFramesNew(i); %#ok<NASGU>
        save( fullfile(dataDir,[ fileName '_' num2str(nFramesNew(i)) 'frames']), 'video', 'fps', 'cmppxx' )
        video = videoOld;
        fps = fpsOld;
    end
end
pwv = zeros(1,length( nFramesNew ));
maxPwv = zeros(1,length( nFramesNew ));
minFps = 1;
maxFps = length( nFramesNew );
% i = 0;
opts.bCorrWaves = 0;
for iPWV = minFps : maxFps
%     i = i+1;
    disp(iPWV)
    [ pwv(iPWV), ~, maxPwv(iPWV) ] = ... 
        PWV( [ fileName '_' num2str(nFramesNew(iPWV)) 'frames'], nPoints, opts );
%     if i>1 
%         plot( [maxPwv(iPWV-1) maxPwv(iPWV)], [pwv(iPWV-1) pwv(iPWV) ],'b-')
%     end
%     hold on
%     drawnow    
end