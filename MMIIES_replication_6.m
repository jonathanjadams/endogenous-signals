

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Macroeconomic Models with Incomplete Information and Endogenous Signals
%  
%                           Jonathan J. Adams
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This version: 6
% Date: 12/19/2024
% Contact: adamsjonathanj@gmail.com

% Dependencies: SIGNAL_OP Package of MATLAB functions, Version 1.5+ (add to path)
% MMIIES_beauty_contest, MMIIES_singleton, MMIIES_confoundingdynamics
% downloadable from jonathanjadams.com



clear all
close all

version = 6;

plot_saving = 0; %switch to 1 to save plots in a subfolder "graphs"


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Uniqueness/Instability in the Confounding Dynamics Models
%                            (Figures 1 and 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MMIIES_confoundingdynamics


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Indeterminancy in the Beauty Contest Model
%                                (Figure 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional dependencies: 
% beauty_contest_discriminant.m 

MMIIES_beauty_contest

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Solutions to the Singleton Model
%                           (Figures 3 and 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional dependencies:
% Matrix of results from the Nimark algorithm: 'nimark_pirfs.mat'
% Matrices of results from Han et al (2022) algorithm: 'ztran_singleton.mat', 'ztran_singleton_time.mat'

MMIIES_singleton



