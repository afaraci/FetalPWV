function [ pwv, maxPwv, aortaL ] = increasingLength( fileName, nPoints )

%%
mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
dataDir = fullfile( mainDir, 'Data' );

%%
fps = [];
cmppxx = [];
load( fullfile( fileName ), 'video', 'fps', 'cmppxx' )

opts.bIncreasingLength = 1;

opts.aortaLength = [];
opts.bCorrWaves = 0;
opts.bDebug = 0;
opts.frameStep = 10;
opts.bDebugWaves = 0;

[ ~, ~, ~, maxPwv, d, dis ] = PWV( fileName, nPoints, opts );

figure
plot( [ 0 maxPwv ], [ -3 -3 ],'Color',[1 0 0],'LineStyle',':','LineWidth', 1 )
hold on
plot( [ 0 maxPwv ], [  0 -maxPwv ], 'Color', [1 0 0], 'LineWidth', 1 )
ax = gca;
axis square
axis tight
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.YLabel.String = 'PWV [m/s]';
ax.XLabel.String = 'Max PWV [m/s]';
ax.FontSize = 12;
ax.FontWeight = 'Bold';
ax.FontName = 'Calibri';
hold on
line( [3 3], ax.YLim, 'LineWidth', 1 )
title(sprintf('Estimated PWV\n vs length of aorta'))

tRes = 1;
pwv = zeros(1,length(d)-1);
maxPwv = pwv;

for i = 2 : length(d)
    i
    [ pwvEst, maxPwv(i-1) ] = ...
        crossCorr( d(1:i), dis(:,1:i), cmppxx, fps, tRes, corrType );
    pwv(i-1) = pwvEst.mean;
    if i > 2
        plot( [ maxPwv(i-2) maxPwv(i-1) ], [ pwv(i-2) pwv(i-1) ], ...
            'Color', [0 0 1], 'Marker', 'o', 'MarkerFaceColor', [1 1 0]);
        drawnow
    end
end
aortaL = maxPwv/fps*1000;
end