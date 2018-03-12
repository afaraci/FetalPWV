function [fitresult, gof, output] = createFitfourier1(t, sig)
%CREATEFIT(T,SIG)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : t
%      Y Output: sig
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 27-Oct-2017 18:42:11


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, sig );

% Set up fittype and options.
ft = fittype( 'fourier1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0 0 15.7612106010607];

% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft, opts );

fitresult.a1=0;
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'sig vs. t', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel t
ylabel sig
grid on

