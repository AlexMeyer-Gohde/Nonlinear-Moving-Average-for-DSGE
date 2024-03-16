function [answer]=own_A_times_B_kronecker_C(A,B,C,options_,MA_,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% own_A_times_B_kronecker_C.m
%
% This file is a wrapper for calculating A*kron(B,C). Currently, the mex
% files from Dynare are used (version detection corrects for the change in
% the order of output in version 4.3.0).
% 
%Included in the package nonlinear_MA
%
%THIS VERSION: 1.0.8 Jun. 3, 2013
%
%Copyright: Hong Lan and Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==5
    if MA_.dynare_version<4.3 %Starting in Dynare 4.3.0, dr1 is removed. 
        [waste,answer] =A_times_B_kronecker_C(A,B,C,options_.threads.kronecker.A_times_B_kronecker_C);
    elseif MA_.dynare_version>=4.3 %replace with stochastic_solvers in 4.3
        [answer, waste] =A_times_B_kronecker_C(A,B,C,options_.threads.kronecker.A_times_B_kronecker_C);
    else
        disp('Error, no certainty evuivalent dr calculated');
        answer=[];
    end
else
         disp('Error, too many arguements passed to own_A_times_B_kronecker_C');
        answer=[];
 end