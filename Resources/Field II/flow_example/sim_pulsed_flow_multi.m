%  Simulate a pulsatile flow from the femoral artery for a spectral velocity system.
%
%  The code can be run on multiple computers in parallel and does not
%  store the scatterers on disk. It can also be restarted if some
%  simulations are not completed. This is done by deleting the
%  affected file.
%
%  This version uses the new functions for set-up
%
%  Version 1.0 by Joergen Arendt Jensen, June 11, 2012.

%  Make a short random pause to not make evryone start at the same time

pause(abs(randn(1)))

%  Set the angle between the ultrasound beam and the velocity vector for this simulation

sim_angle=70;

%  This commands should start the Field II simulation program

field

%  Setup the standard simulation parameters and make the transducer

[sys, focus_para] = set_standard_parameters();
[xdc] = define_convex_transducer (sys, focus_para);

%  Set overall parameters for the simulation

%  Physical parameters for the femoral artery

sim.rho=1.06e3;          %  Density of blood   [kg/m^3]
sim.mu=0.004;            %  Viscosity of blood [kg/m s]
sim.R=4.2/1000;          %  Radius of vessel   [m]
sim.omega_0=2*pi*60/60;  %  Angular frequency, heart rate = 60 beats/min
sim.V0=0.15;             %  Mean velocity      [m/s]
sim.Nt=100;              %  Number of time values
sim.deltaR=0.002;        %  Radial sampling interval
sim.r_rel=(0:sim.deltaR:1)';  %  Values for the relative radius

%   Calculate psi

calc_psi

%   Data for phase and amplitude of waveform

sim.Vp=[1.00    0.00
    1.89   32
    2.49   85
    1.28  156
    0.32  193
    0.27  133
    0.32  155
    0.28  195
    0.01  310];

%  Convert the data to the correct amplitude and to radians

sim.Vp(:,1)=sim.Vp(:,1)*sim.V0;
sim.Vp(:,2)=sim.Vp(:,2)*pi/180;

%  Set parameters for flow simulation

sim.Tprf=1/15000;             %  Pulse repetition frequency  [s]
sim.N_shoots=ceil(1/sim.Tprf);     %  Number of shoots for one second

sim.x_range=40/1000;         %  Range of scatterers in x-range [m]
sim.y_range=2.5*sim.R;       %  Range of scatterers in y-range [m]
sim.z_range=2.5*sim.R;       %  Range of scatterers in z-range [m]
sim.z_offset=40/1000-sim.R;  %  Distance of center of vessel from transducer [m]
sim.depth=40/1000;           %  Center of vessel and place for optimal TO focusing

%  Calculate the number of scatterers needed from the estimated size of
%  the PSF. These should be found from the point spread fucntion.
%  There should be 10 scatterers within a resolution cell in order to
%  get fully developed speckle

fwhm_ax=xdc.lambda;      %  Axial resolution [m]
fwhm_lat=2*xdc.lambda;   %  Lateral resolution [m]
fwhm_az=5*xdc.lambda;    %  Azimuthal resolution [m]

sim.N_scat=round(10*sim.x_range*sim.y_range*sim.z_range/ ...
    (fwhm_lat   *fwhm_az    *fwhm_ax)     );             %  Number of scatterers in phantom
disp(['Number of scatterers in the phantom calculated to ',num2str(sim.N_scat)])

%  Determine directory name and make it, if it does not exits

dir_name='rf_data/';
if (exist(dir_name) ~= 7)  %#ok<*EXIST> %  Create directory
    mkdir(dir_name)
end

%  Make for the simulation

sim.theta=(90-sim_angle)/180*pi;         %  Angle between ultrasound beam and velocity vector.

%  Directory name for storage

dir=[dir_name,'rf_data_',num2str(sim_angle),'_deg/'];  %  Directory and name for storage
if (exist(dir) ~= 7)  %  Create directory
    mkdir(dir)
end

%  Find the name of this machine

[dummy, hostname] = system('hostname');

%  Set the display time to 20 seconds

set_field('show_times',20);

%  Set the seed of the generator, so we get the
%  same "random" data

% reset(RandStream.getDefaultStream);
reset(RandStream.getGlobalStream); % Changed by AF 03/04/17

%  Generate the coordinates and amplitude
%  Coordinates are rectangular within the range.
%  The amplitude has a Gaussian distribution.

x=sim.x_range*(rand(1,sim.N_scat)-0.5);
y=sim.y_range*(rand(1,sim.N_scat)-0.5);
z=sim.z_range*(rand(1,sim.N_scat)-0.5);

%  Find which scatterers that lie within the blood vessel

r=(y.^2+z.^2).^0.5;
within_vessel= (r < sim.R);

%  Assign an amplitude for each scatterer
%  Set the amplitude outside of the vessel to 10 times the inside

amp=randn(sim.N_scat,1).*within_vessel'+10*randn(sim.N_scat,1).*(1-within_vessel');

%  Do for a all the shoots

for k=1:sim.N_shoots
    
    %  Find the time
    
    time=k*sim.Tprf;
    
    %  Calculate the profile
    
    disp (['Time is ',num2str(time),' seconds'])
    prof=2*sim.V0*(1-sim.r_rel.^2);
    index = 1 : max(size(sim.r_rel));
    for p=2:8
        prof=prof + sim.Vp(p,1)*abs(psi(index,p-1)).*cos((p-1)*sim.omega_0*time - sim.Vp(p,2) + angle(psi(index,p-1)));
    end;
    
    %  Generate the rotated and offset block of scatterers
    
    xnew=x*cos(sim.theta)+z*sin(sim.theta);
    znew=z*cos(sim.theta)-x*sin(sim.theta) + sim.z_offset;
    positions=[xnew; y; znew;]';
    
    %   Calculate the received signals for all elements
    %   if it has not been done before
    
    %  Check if the file already exits
    
    file_name = [dir,'rf_ln_',num2str(k),'.mat'];
    cmd = 'fid = fopen(file_name,''r'');';
    eval(cmd);
    
    %  Do the processing if the file does not exits
    
    if (fid == -1)
        cmd=['save ',dir,'rf_ln_',num2str(k),'.mat k hostname'];
        eval(cmd);
        
        %  Save the parameters for the simulation
        
        if (k==1)
            
            cmd=['save ',dir,'parameters.mat sim sys xdc focus_para'];
            eval(cmd);
        end
        
        %   Calculate the received response
        
        [rf_data, t_start]=calc_scat(xdc.transmit_aperture, xdc.receive_aperture, positions, amp);
        
        %  Store the result
        
        cmd=['save ',dir,'rf_ln_',num2str(k),'.mat rf_data t_start'];
        eval(cmd)
    else
        fclose (fid);
    end
    
    %  Find the velocity of the scatterers for this time and position
    
    in_velo = floor(r/(sim.R*sim.deltaR).*within_vessel) + 1;
    velocity = prof(in_velo)'.*within_vessel;
    
    %  Propagate the scatterers and aliaze them
    %  to lie within the correct range
    
    x=x + velocity*sim.Tprf;
    outside_range= (x > sim.x_range/2);
    x=x - sim.x_range*outside_range;
    outside_range= (x < -sim.x_range/2);
    x=x + sim.x_range*outside_range;
end