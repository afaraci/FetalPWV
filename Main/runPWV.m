function [ pwvX, pwvP, pwvF, d ] = runPWV( fileName, nPoints, opts )

if nargin < 3, opts = []; end
if isfield(opts,'bDebug'),          opts.bDebug          = opts.bDebug;          else, opts.bDebug          = 1; end
if isfield(opts,'bDebugWaves'),     opts.bDebugWaves     = opts.bDebugWaves;     else, opts.bDebugWaves     = 1; end
if isfield(opts,'bPixelSelection'), opts.bPixelSelection = opts.bPixelSelection; else, opts.bPixelSelection = 1; end
if isfield(opts,'bAskRot'),         opts.bAskRot         = opts.bAskRot;         else, opts.bAskRot         = 1; end
if isfield(opts,'bRot'),            opts.bRot            = opts.bRot;            else, opts.bRot            = 1; end
if isfield(opts,'bCorrWaves'),      opts.bCorrWaves      = opts.bCorrWaves;      else, opts.bCorrWaves      = 1; end
if isfield(opts,'bPlotBoxes'),      opts.bPlotBoxes      = opts.bPlotBoxes;      else, opts.bPlotBoxes      = 0; end

close all
[ pwvX, pwvP, pwvF, d ] = PWV( fileName, nPoints, opts );