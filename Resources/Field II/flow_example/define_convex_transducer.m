%  Function to set-up the transducer for simulation of a convex array 
%
%  Field II initialization should be done before calling
%
%  The aperture are defined. Note that no initial focusing is set-up
%
%  Calling:  [xdc] = define_convex_transducer (sys, focus_para); 
% 
%   Parameters:  sys        -  Structure with system parameters
%                focus_para - Structure definint the focusing parameters
% 
%   Return:      xdc  -  Structure with all definitions for the transducer
% 
%  Version 1.1 by Joergen Arendt Jensen, June 11, 2012.

function [xdc] = define_convex_transducer (sys, focus_para) 

%  Set the initial parameters for the transducer
%  The structure xdc contains all parameters for the transducer

xdc.f0=3.42e6;                  %  Transducer center frequency [Hz]
xdc.M0=4;                       %  Number of cycles in emitted pulse

xdc.lambda=sys.c/xdc.f0;        %  Wavelength [m]

xdc.pitch=0.33/1000;            %  Element pitch [m]
xdc.kerf=xdc.pitch/20;          %  Kerf - actually not indicated in the data sheet [m]
xdc.width=xdc.pitch-xdc.kerf;   %  Width of element
xdc.height=13/1000;             %  Height of element [m]
xdc.elevation_focus=90/1000;    %  Fixed elevation focus [m]
xdc.N_elements=192;             %  Number of physical elements of the transducer
xdc.Rconvex=60/1000;            %  Convex radius [m]

xdc.focus_transmit=focus_para.focus_transmit;              %  Electronic focus in transmit [m]
xdc.F_number_transmit=focus_para.F_number_transmit;        %  F# in transmit           
xdc.N_e_receive=focus_para.N_e_receive;                    %  Number of receive elements

xdc.depth=focus_para.depth;                                %  Depth for optimizing the receive focusing [m]

% Define the transducer i both transmit and receive

xdc.focus_transmit_point=[0 0 xdc.focus_transmit];      %  Fixed emit focal point [m]
xdc.transmit_aperture = xdc_convex_focused_array (xdc.N_elements, xdc.width, xdc.height, xdc.kerf, ...
                            xdc.Rconvex, xdc.elevation_focus, 2, 10, xdc.focus_transmit_point);
xdc.receive_aperture  = xdc_convex_focused_array (xdc.N_elements, xdc.width, xdc.height, xdc.kerf, ...
                            xdc.Rconvex, xdc.elevation_focus, 2, 10, [0 0 xdc.depth]);

%  Set the impulse response and excitation of the aperture

impulse_response=sin(2*pi*xdc.f0*(0:1/sys.fs:2.5/xdc.f0));
xdc.impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (xdc.transmit_aperture, xdc.impulse_response);

xdc.excitation=sin(2*pi*xdc.f0*(0:1/sys.fs:xdc.M0/xdc.f0));
xdc.excitation=xdc.excitation .* hanning(length(xdc.excitation))';
xdc_excitation (xdc.transmit_aperture, xdc.excitation);

%  Set the impulse response for the receive apertures

xdc_impulse (xdc.receive_aperture, xdc.impulse_response);

%  Calculate the number of elements in transmit from the F# in transmit

aperture_width=xdc.depth/xdc.F_number_transmit;        %  Optimized aperture size for transmit
xdc.N_e_transmit=round(aperture_width/xdc.pitch);      %  Number of transmit elements
if rem(xdc.N_e_transmit,2)==1                          %  Ensure that an odd number of elements is used during transmit
   xdc.N_e_transmit=xdc.N_e_transmit+1;
   end
xdc.apo_receive=[zeros(1, round((xdc.N_elements-xdc.N_e_receive)/2)) hamming(xdc.N_e_receive)' zeros(1, round((xdc.N_elements-xdc.N_e_receive)/2))];

apo2=hamming(xdc.N_e_transmit)';
xdc.apo_transmit=[zeros(1, round((xdc.N_elements-xdc.N_e_transmit)/2)) apo2 zeros(1, round((xdc.N_elements-xdc.N_e_transmit)/2))];

% Define the transducer in both transmit and receive

xdc.focus_transmit_point=[0 0 xdc.focus_transmit];      %  Fixed emit focal point [m]
xdc_focus(xdc.transmit_aperture, 0, xdc.focus_transmit_point);

%  Set the receive apodization function so that a transverse oscillation is made

xdc_apodization (xdc.receive_aperture, 0, xdc.apo_receive);

%  Only use the center elements for transmit

xdc_apodization (xdc.transmit_aperture, 0, xdc.apo_transmit);