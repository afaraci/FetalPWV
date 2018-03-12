function [ Pt, phantom_amplitudes ] = generatePWVphantom( pwv, fps, N, bDebug )

if nargin < 1;    pwv =   -3; end
if nargin < 2;    fps =  148; end % Frames per second
if nargin < 3;      N = 1000; end % Number of scatterers in phantom
if nargin < 4; bDebug =    0; end

mainDir = cd;
mainDir = mainDir(1:strfind(mainDir,'FetalPWV')+length('FetalPWV')-1);
dataDir = fullfile( mainDir, 'Data' );

load( fullfile(dataDir,'simsig'),'fittedmodel2')
period = 2*pi/fittedmodel2.w; % wave period [s]
nFrames = round( period*fps );
t = 0 : period/(nFrames-1) : period;
t(end) = [];
nFrames = nFrames - 1;
%  Load the computer phantom data
[ phantom_positions, phantom_amplitudes ] = tissue_pht( N, 1 );

phantom_positions(:,2) = 0; % For 2D phantom

%%  Phantom dimensions
x_size = 40/1000;            % Width of phantom [m]
% y_size = y_size/1000;        % Transverse width of phantom [m]
z_size = 20/1000;            % Height of phantom [m]
z_plus = z_size/2 + 50/1000; % Start of phantom surface [m];

%%
% lenAor = x_size*1000; % [mm]
R = fittedmodel2.a0 / 2 / 1000; %  Radius of blood vessel [m]
Wt = R/10; % Thickness of vessel's wall
y = phantom_positions(:,2);
z = phantom_positions(:,3);
% The inside includes the wall
inside = y.^2 + (z-z_plus).^2 <= (R + Wt)^2;
outside = ~inside;
phantom_positions_out = phantom_positions(outside,:);
phantom_amplitudes = phantom_amplitudes(outside);

%% Vessel's wall
df = 30; % density factor
% Load the computer phantom data
% Concentration of points in the vessel's wall is df-times higher than in
% the rest of the phantom's tissue
[ phantom_positions2, phantom_amplitudes2 ] = tissue_pht( N*df, 1 );
y = phantom_positions2(:,2);
z = phantom_positions2(:,3);
wall = logical( ( y.^2 + (z-z_plus).^2 >=        R^2 ) .* ... 
                ( y.^2 + (z-z_plus).^2 <= (R + Wt)^2 )    );
% figure
% plot3(phantom_positions_out(:,1),phantom_positions_out(:,2),phantom_positions_out(:,3),'b.')
% hold on
% plot3(phantom_positions2(wall,1),phantom_positions2(wall,2),phantom_positions2(wall,3),'b.')
% campos( [ -0.0000   -0.2235    0.0600 ] )
% axis equal
% axis tight

phantom_positions_wall = phantom_positions2(wall,:);
phantom_positions_out = [ phantom_positions_out; phantom_positions_wall ];
phantom_amplitudes2 = phantom_amplitudes2(wall);
phantom_amplitudes  = [ phantom_amplitudes; phantom_amplitudes2 ];
% hold on
% plot3(phantom_positions_out(:,1),phantom_positions_out(:,2),phantom_positions_out(:,3),'b.')

%%
for iP = 1 : length( phantom_positions_out )
    P = phantom_positions_out(iP,:);
    % Radius over time at position P(1)
    Rt = feval( fittedmodel2, t + P(1)/pwv ) / 2;
    sF = 1 - ( abs( P(3) - z_plus ) - R ) / ( z_size/2 - R );
    zt = P(3)*1000 + (Rt-R*1000)*sF * ( P(3) - z_plus )/abs( P(3) - z_plus );
    for iFrame = 1 : nFrames
        Pt(iP,:,iFrame) = [ P(1) P(2) zt(iFrame)/1000 ];
    end
end

Pt(:,2,:) = 0; % For 2D phantom

%%
if bDebug
    close all
    figure
    for iFrame = 1 : nFrames
        plot3(Pt(:,1,iFrame),Pt(:,2,iFrame),Pt(:,3,iFrame),'r.')
        campos( [ -0.0000   -0.2235    0.0600 ] )
        axis equal
        axis tight
        drawnow
    end
end