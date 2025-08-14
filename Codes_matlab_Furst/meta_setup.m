% Matlab Script to assemble Favier model results. Start plotter.m afterwards
%
% script by Johannes J. Fuerst, 2013
%
% INPUT:        see scripts in ./bin/
%               read_PIG_favier.m
%               input_pig.mat (produce with ../Input_bm2_db_vel/meta_region.m)
%               ../../../input/favier2014/*/*/export*.csv
%
% OUTPUT:       figures
%

% by Johannes J. Fuerst

clear all
close all

path(path,'./bin/')

% =====================================================================
% CONSTANTS
% =====================================================================

grav     = 9.81;
rhoice   = 917;
rhowater = 1028;
sealevel = 0;
MPa      = 1e6;
secy     = 365.25*24*60*60;

nexp     = 24;

% =====================================================================
% VARIABLES
% =====================================================================

% =====================================================================
% LOAD ALL 
% =====================================================================

% get bedmap2 information on PIG
tic
display('BEDMP2')
load('../Input_bm2_db_vel/data/input_pig.mat');
toc

% load PIG - only surface data
tic
display('Calving experiment data - Favier et al. (2014)')
read_PIG_favier
toc

% calculate buttressing on Elmer nodes from surface fields
tic
display('Determine buttressing on Elmer nodes')
buttressing
toc

% interpolate buttressing and velocities on BEDMAP2
tic
display('Interpolate to BEDMAP2 grid')
regridder
toc

% =====================================================================
% STORAGE 
% =====================================================================

%save('./data/insar.mat','insar','-v7.3')
%save('./data/bedmap2.mat','bedmap2','-v7.3')
%clear insar
%}












