%  Make a the scatterers for a simulation of the rf signal from
%  flow in the femoral artery.
%  The result is stored as scatterer positions in a file
%
%  Version 1.0 9/6-97, JAJ

%  Set the seed of the random number generator

randn('seed',sum(100*clock))

%  Physical parameters for the femoral artery

c=1540;              %  Ultrasound speed of sound [m/s]
f0=2e6;              %  Center frequency of transducer [Hz] 
rho=1.06e3;          %  Density of blood   [kg/m^3]
mu=0.004;            %  Viscosity of blood [kg/m s]
R=4.2/1000;          %  Radius of vessel   [m]
omega_0=2*pi*62/60;  %  Angular frequency, heart rate = 62 beats/min
V0=0.15;             %  Mean velocity      [m/s]
Nt=100;              %  Number of time values
deltaR=0.002;        %  Radial sampling interval

r_rel=(0:deltaR:1)';  %  Values for the relative radius

%   Calculate psi

calc_psi

%   Data for phase and amplitude of waveform

Vp=[1.00    0.00
    1.89   32 
    2.49   85  
    1.28  156
    0.32  193 
    0.27  133 
    0.32  155
    0.28  195 
    0.01  310];

%  Convert the data to the correct amplitude and to radians

Vp(:,1)=Vp(:,1)*V0;
Vp(:,2)=Vp(:,2)*pi/180;

%  Initialize the ranges for the scatteres
%  Notice that the coordinates are in meters

x_range=0.04;    %  x range for the scatterers  [m]
y_range=2*R;     %  y range for the scatterers  [m]
z_range=2*R;     %  z range for the scatterers  [m]
z_offset=0.06;   %  Offset of the mid-point of the scatterers [m]

%  Set the number of scatterers. It should be roughly
%  10 scatterers per cubic wavelength

lambda=c/f0;
N=floor(10*x_range/(5*lambda)*y_range/(5*lambda)*z_range/(lambda*2));
disp([num2str(N),' Scatterers'])

%  Generate the coordinates and amplitude
%  Coordinates are rectangular within the range.
%  The amplitude has a Gaussian distribution.

x=x_range*(rand(1,N)-0.5);
y=y_range*(rand(1,N)-0.5);
z=z_range*(rand(1,N)-0.5);

%  Find which scatterers that lie within the blood vessel

r=(y.^2+z.^2).^0.5;
within_vessel= r < R;

%  Assign an amplitude for each scatterer

blood_to_stationary= 10;   %  Ratio between amplitude of blood to stationary tissue
amp=randn(1,N).*((1-within_vessel) + within_vessel*blood_to_stationary);
amp=amp';

%  Generate files for the scatteres over a number of pulse emissions

Tprf=1/5e3;   %  Time between pulse emissions  [s]
Nshoots=5000; %  Number of shoots

i=0;
for i=1:Nshoots
i
  time=i*Tprf;

  %  Calculate the profile

  disp (['Time is ',num2str(time),' seconds'])
  prof=2*V0*(1-r_rel.^2);
  index=[1:max(size(r_rel))];
  for p=2:8
    prof=prof + Vp(p,1)*abs(psi(index,p-1)).*cos((p-1)*omega_0*time - Vp(p,2) + angle(psi(index,p-1)));
    end;


  %  Generate the rotated and offset block of sample

  theta=60/180*pi;
  xnew=x*cos(theta)+z*sin(theta);
  znew=z*cos(theta)-x*sin(theta) + z_offset;
  positions=[xnew; y; znew;]';

  %   Save the matrix with the values

  cmd = ['save /data3/users_nbu/jaj/sim_flow_data/scat_data/scat_',num2str(i),'.mat positions amp']
  eval(cmd)

  %  Find the velocity of the scatterers for this time and position

  in_velo = floor(r/(R*deltaR).*within_vessel) + 1;
  velocity = prof(in_velo)'.*within_vessel;
  plot(r,velocity,'x')
  drawnow

  %  Propagate the scatterers and aliaze them
  %  to lie within the correct range
  
  x1=x;
  x=x + velocity*Tprf;
  outside_range= (x > x_range/2);
  x=x - x_range*outside_range;
  end
  


