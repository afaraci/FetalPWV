%  Set standard parameter values for the simulation. 
%
%  Calling: [sys focus_para] = set_standard_parameters(); 
% 
%   Parameters:  None
% 
%   Return:      sys         -  Structure with system parameters
%                focus_para  -  Structure with all definitions for the focusing
% 
%  Version 1.1 by Joergen Arendt Jensen, June 11, 2012.

function [sys focus_para] = set_standard_parameters() 

%  The structure sys contains all parameters for the simulation system

sys.fs=100e6;                   %  Sampling frequency [Hz]
sys.c=1540;                     %  Speed of sound [m/s]

%  Set the sampling frequency and display time

set_sampling(sys.fs);
set_field('c', sys.c);

%  Set different values for experimentation with focusing

focus_para.depth=40/1000;            %  Depth for optimizing the receive focusing [m]
focus_para.focus_transmit=50/1000;   %  Electronic focus in transmit [m]
focus_para.F_number_transmit=2;      %  F# in transmit (used for setting the number of elements)          
focus_para.N_e_receive=64;           %  Number of receive elements