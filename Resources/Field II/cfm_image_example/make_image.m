function [ new_env, xSpace, ySpace ] = ...
    make_image( rf_lines, tstart_lines, image_width )

%  Compress the data to show 40 dB of
%  dynamic range for the CFM phantom image
%
%  version 1.2 by Joergen Arendt Jensen, April 6, 1997.

%% Parameters
f0 = 3.5e6; %  Transducer center frequency [Hz]
fs = 100e6; %  Sampling frequency [Hz]
c = 1540;  %  Speed of sound [m/s]

%% Load RF data lines and adjust them in time
min_sample = 0;
nLines = length(rf_lines);      % Number of lines in image
image_width = image_width/1000; % Size of image sector [m]
d_x = image_width/nLines;       % Increment for image  [m/line]

%% Find envelope
rf_data = cell2mat( rf_lines( 1 ) );
tstart  = tstart_lines( 1 );
envLength = tstart*fs + length(rf_data);
env = zeros( envLength, nLines );
parfor iLine = 1 : nLines
    rf_data = cell2mat( rf_lines( iLine ) );
    tstart  = tstart_lines( iLine );
    if tstart > 0
        rf_env = abs( hilbert( [ zeros( round(tstart*fs-min_sample),1); rf_data ] ));
    elseif tstart == 0
        rf_env = abs( hilbert( rf_data( abs(tstart*fs)+1 : length(rf_data) ) ));
    end
    rf_env(end+1 : envLength) = zeros;
    env( :, iLine ) = rf_env;
    % env( 1 : length(rf_env), iLine ) = rf_env;
    % env( iLine ) = { rf_env };
end

%%  Logarithmic compression
env = env - min( env(:) );
env = env / max( env(:) );
new_env = log( env + 0.01 );
D = 1;        % Sampling frequency decimation factor
ID_bmode = 1;
% new_env = imresize( new_env, [ floor(n/D), m*ID_bmode ] );
new_env = new_env - min( new_env(:) );
new_env = new_env / max( new_env(:) );
n = size( new_env, 1 );
fn_bmode = fs/D;

%% Display B-mode image
xSpace = ( ( 1:(ID_bmode*nLines-1) )*d_x/ID_bmode-nLines*d_x/2 )*1000; % [mm]
ySpace = ( (1:n)/fn_bmode + min_sample/fs )*c/2*1000; % [mm]

% image( xSpace, ySpace, new_env )
% xlabel('Lateral distance [mm]')
% ylabel('Axial distance [mm]')
% colormap( gray(nColours) )
% axis image
% axis equal
% xlim = image_width*1000/2; % [mm]
% axis( [ -xlim xlim 51 71 ] )