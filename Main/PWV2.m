function [ pwvP, pwvF ] = PWV2( w, fps, cmppxx, d )

w = w * cmppxx * 10; % [mm]
nF = size( w, 1 );
nW = size( w, 2 );

t = ( 1 : nF ) / fps; % [s]
z1 = cell(1,nW);
for iW = 1 : nW
    sp(iW) = spaps( t, w(:,iW), 0, 3 );
    fnz1 = fnzeros( fnder( sp(iW) ) );
    fnz1(2,:) = [];
    z1( iW ) = { fnz1 };
end

for i = 1 : length(z1)-1
    [~, idx ] = max( xcorr(cell2mat(z1(i)), cell2mat(z1(i+1))) );
    a = cell2mat(z1(i));
    b = cell2mat(z1(i+1));
    z1(i)   = { a(end-idx+1:end) };
    z1(i+1) = { b(end-idx+1:end) };
end
z1 = reshape(cell2mat(z1),length(cell2mat(z1))/length(z1),length(z1))';

val1 = zeros( nW, size(z1,2) );
for iW = 1 : nW
    val1(iW,:) = fnval( z1(iW,:), sp(iW) );
    fnz2 = fnzeros( fnder( sp(iW), 2 ), [z1(iW,1) z1(iW,end)] );
    z2(iW,:) = fnz2(1,:);
    val2(iW,:) = fnval( z2(iW,:), fnder(sp(iW)) );
end

ttp = z1(:,val1(1,:) > 0); % Time-to-peak
ttb = z1(:,val1(1,:) < 0); % Time-to-bottom
z2 = z2(:,val2(1,:) > 0);
% val2( :, val2(1,:) < 0 ) = [];
iP = zeros( nW, size(z2,2) );
ttf = ones( size( ttb ) );
for iW = 1 : nW
    a = fnval( z2(iW,:), fnder( sp(iW) ) );
    % Inflection points, points of maximum velocity
    iP(iW,:) = fnval( z2(iW,:), sp(iW) ); 
    b = iP(iW,:) - a.*z2(iW,:);
    % y = fnval( ttb(iW,:), sp(iW) );
    y = 0;
    ttf(iW,:) = ( y - b )./a;
%     fnz3 = fnzeros( fnder( sp(iW), 3 ), [ttf1(iW,1) ttp(iW,end)] );
%     ttf3(iW,:) = { fnz3(1,:) };
end

% for i = 1 : length(ttf3)-1
%     [~, idx ] = max( xcorr(cell2mat(ttf3(i)), cell2mat(ttf3(i+1))) );
%     a = cell2mat(ttf3(i));
%     b = cell2mat(ttf3(i+1));
%     ttf3(i)   = { a(end-idx+1:end) };
%     ttf3(i+1) = { b(end-idx+1:end) };
% end
% ttf3 = reshape(cell2mat(ttf3),length(cell2mat(ttf3))/length(ttf3),length(ttf3))';

dm = d*cmppxx*10/1000; % [m]

pP  = polyfit( dm, mean( ttp-ttp(1,:), 2 ), 1 );
pF  = polyfit( dm, mean( ttf-ttf(1,:), 2 ), 1 );
% pF2 = polyfit( dm, mean( ttf1-ttf1(1,:), 2 ), 1 );
% pF3 = polyfit( dm, mean( ttf3-ttf3(1,:), 2 ), 1 );

pwvP  = ( 1/pP(1) ); % [m/s]
pwvF  = ( 1/pF(1) ); % [m/s]
% pwvF2 = ( 1/pF2(1) ); % [m/s]
% pwvF3 = ( 1/pF3(1) ); % [m/s]