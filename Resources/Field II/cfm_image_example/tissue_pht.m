%  Create a computer model of the tissue surrounding a
%  vessel
%
%  Calling: [positions, amp] = tissue_pht (N);
%
%  Parameters:  N - Number of scatterers in the phantom
%
%  Output:      positions  - Positions of the scatterers.
%               amp        - Amplitude of the scatterers.
%
%  Version 1.0, March 26, 1997 by Joergen Arendt Jensen

function  [ positions, amp ] = tissue_pht( N, y_size )

if nargin < 1
    N = 1000;
end
if nargin < 2
    y_size = 1; % [mm]
end

%%  Phantom dimensions
x_size = 40/1000;            % Width of phantom [m]
y_size = y_size/1000;        % Transverse width of phantom [m]
z_size = 20/1000;            % Height of phantom [m]
z_plus = z_size/2 + 50/1000; % Start of phantom surface [m];

%%  Initialise vessel ranges
R = 2/1000;         %  Radius of blood vessel [m]
y_range = 2*R;      %  y range for the scatterers  [m]
z_range = 2*R;      %  z range for the scatterers  [m]
z_offset = 60/1000; %  Offset of the mid-point of the scatterers [m]

%%  Create the general scatterers
x = ( rand(N,1)-0.5 )*x_size;
y = ( rand(N,1)-0.5 )*y_size;
z = ( rand(N,1)-0.5 )*z_size + z_plus - z_offset;

%%  Generate the amplitudes with a Gaussian distribution
amp = randn(N,1);

%%  Make the vessel
inside = y.^2 + z.^2 < R^2;
outside = 1-inside;
% Reduce inside amplitudes by 40 dB
amp = amp.*outside + amp.*inside/db2mag(40);

%%  Generate the offset block of sample
z = z + z_offset;

%%  Return the variables
positions = [ x y z ];