%% Import Script for ODF Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = crystalSymmetry('m-3m', [1 1 1], 'color', 'light blue');

% specimen symmetry
SS = specimenSymmetry('1');

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

% path to files
pname = 'D:\data\ODF';

% which files to be imported
fname = [pname '\test.txto'];

%% Import the Data

% specify kernel
psi = deLaValeePoussinKernel('halfwidth',10*degree);

% create an EBSD variable containing the data
odf = loadODF(fname,CS,SS,'density','kernel',psi,'resolution',5*degree,...
  'interface','generic',...
  'ColumnNames', { 'phi1' 'Phi' 'phi2'}, 'Columns', [1 2 3], 'Bunge');

plotPDF(odf,[Miller(0,0,0,1,CS),Miller(0,1,-1,1,CS),Miller(1,1,-2,1,CS)])

CLim(gcm,[0,2.5])

mtexColorbar