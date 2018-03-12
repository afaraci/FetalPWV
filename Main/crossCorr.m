function [ pwv, maxPwv ] = crossCorr( fileName, d, dis, cmppxx, fps, ...
    bIncreasingLength, bDebug )

if nargin < 6, bIncreasingLength = 0; end
if nargin < 7, bDebug = 0; end

%%
mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
if ~( contains(fileName,'simsig') || contains(fileName,'NLdata') ...
        || contains(fileName,'phantom') )
    dataDir = fullfile( mainDir, 'Data', 'Clinical' );
else
    dataDir = fullfile( mainDir, 'Data', 'Phantom' );
end

%% Initilisation
% e = 3 : 1022;
e = d(1) : d(end);
e(d-d(1)+1) = [];
maxPwv = (d(end)-d(1))*cmppxx*10*fps/1000;
fprintf( '\nMax PWV = %1.2f m/s\n', maxPwv );

%%
load( fullfile( dataDir, fileName ), 'mm', 'd1' );

%% Computing the cross correlation matrix
% cond = ~exist('mm','var');
cond = 1;
if cond
    mm  =  zeros( size(dis,2), size(dis,2) );
    for j = 1 : size( dis, 2 )
        m = dis(:,j:end);
        mm( j, j:size( dis, 2) ) = corr( m(2:end,1),m(1:end-1,:) );
    end
    d1 = d(1); %#ok<NASGU>
    save( fullfile( dataDir, fileName ), 'd1', 'mm', '-append' );
end

% mm(:,e-2) = [];
% mm(e-2,:) = [];

numLengths =  1;
if bIncreasingLength, numLengths = length(d)-1; end
pwv = zeros(1,numLengths);

%% Interpolating cross correlation
for iP = 1 : numLengths
    if numLengths > 1, disp( iP ), end
    di = d(1):d(end);
    mmi = zeros( size(mm,1), length(di) );
    
    for j = 1 : size(mmi,1)-1
        dT = d(j:end);
        mT = mm(j,j:end);
        dTi = d(j):d(end);
        mTi = interp1( dT, mT, dTi, 'spline' );
        mmi(j,1:length(mTi)) = mTi;
    end
    
    %%
    maxAcori  = zeros(1,length(mmi));
    for j = 1 : size(mmi,2)
        maxAcori(j)  = mean( nonzeros( mmi(:,j)  ) );
    end
    % save( fullfile( dataDir, fileName ), 'maxAcori', '-append' );
    %%
    if bDebug
        figure
        plot( di, maxAcori, '*' );
    end
    
    %%
    cond = 1;
    span = 1;
    while cond
        maxAcoriS = smooth(maxAcori,span,'loess');
        [pks,locs] = findpeaks( maxAcoriS, di );
        cond = sum(pks>0.85)>1;
        span = span+2;
    end
    hold on
    findpeaks( maxAcoriS, di );
%     [~, LOCmax] = max( maxAcori );
    % [~, LOCmin] = min( maxAcori(1:LOCmax) );
%     LOCmin = 1;
%     dd = [ di(LOCmin) di(LOCmax) ];
    dd = [ di(1) locs(pks>0.85) ];
    tt = [ 0 -1 ];
    p = polyfit( dd*cmppxx*10, tt/fps, 1 ); % [s/mm]
    pwv(iP) = ( 1/p(1) )/1000; % [m/s]
    d(end) = [];
    mm(:,end) = [];
    mm(end,:) = [];
end

pwv = flip(pwv);