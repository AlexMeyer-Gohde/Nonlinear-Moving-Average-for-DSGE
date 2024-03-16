function [MA_]=nonlinear_MA(M_,oo_,options_,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nonlinear_MA.m
%
% This is the main program of the package nonlinear_MA
% 
% Note that the folder containing the package must be included in your path
%
% To run, first run dynare with a .mod file and then enter at the command 
% prompt or place at the end of .mod file the following command
%
%  nonlinear_MA(M_,oo_,options_)
%
% The inputs
% M_
% oo_
% options_
% are produced by dynare.
%
% If you would like to change some of the defaults, 
% create a structured array with fields named commesurately with the names
% found below in the section "setup parameters". Pass the structured array
% as a fourth input when calling the program:
%
% nonlinear_MA(M_,oo_,options_, name)
%
% where "name" is whatever name you give to your structured array. For
% example, to change the default of a postive one standard degree shock for
% impulse response to minus three standard degree shock, type 
% >> myoptions.shock_scale=-3;
% and then call nonlinear_MA
% >> nonlinear_MA(M_,oo_,options_, myoptions)
%
% If you would like to pass the shortllist of variables to be plotted
% Dynare's var_list_, simply pass the variable var_list_ as a fifth imput
%
% nonlinear_MA(M_,oo_,options_, [], var_list_)
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
toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 5 %Improper number of inputs
        disp('Too many inputs passed to nonlinear_MA')
else %Otherwise correct number of inputs.
    %I'll start by setting the defaults. If a structure array was passed, 
    %I'll read for override values afterwards
    MA_.shock_scale=1; %Positive one standard deviation shock in IRFs
     MA_.plot_kernels=0;MA_.calculate_kernels=0;MA_.calculate_moments=0;
    if options_.periods>0
        MA_.calculate_simulations=1;
        MA_.plot_simulations=1;
    else
        MA_.calculate_simulations=0;
        MA_.plot_simulations=0;
    end
    if options_.irf>0
        MA_.calculate_irf=1;
        MA_.plot_irf=1;
    else
        MA_.calculate_irf=0;
        MA_.plot_irf=0;
    end
    if nargin>=4 %This is where I'll start overriding the defaults, given
            %that the user passed a strutural array to the program
            %containing override values
        MA_options = varargin{1};
        if isfield(MA_options, 'shock_scale')==1
            MA_.shock_scale=MA_options.shock_scale;
        end
        if isfield(MA_options, 'calculate_kernels')==1
            MA_.calculate_kernels=MA_options.calculate_kernels;
        end
        if isfield(MA_options, 'plot_kernels')==1
            MA_.plot_kernels=MA_options.plot_kernels;
        end
        if isfield(MA_options, 'max_kernel_length')==1
            MA_.max_kernel_length=MA_options.max_kernel_length;
        end
        if isfield(MA_options, 'calculate_simulations')==1
            MA_.calculate_simulations=MA_options.calculate_simulations;
        end
       if isfield(MA_options, 'plot_simulations')==1
            MA_.plot_simulations=MA_options.plot_simulations;
        end
        if isfield(MA_options, 'calculate_irf')==1
            MA_.calculate_irf=MA_options.calculate_irf;
        end
        if isfield(MA_options, 'plot_irf')==1
            MA_.plot_irf=MA_options.plot_irf;
        end
    end
end
if  MA_.calculate_kernels==1
    if  isfield(MA_, 'max_kernel_length')==0
        MA_.max_kernel_length=40;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect model information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract a numeric version of dynare's version
[MA_.dynare_version]=return_dynare_version(dynare_version);
% Definitions from Dynare's dr1.m
MA_.klen = M_.maximum_endo_lag + M_.maximum_endo_lead + 1;
MA_.iyv  = M_.lead_lag_incidence';
MA_.iyv  = MA_.iyv(:);
it_= M_.maximum_lag + 1;

% prepare steady state values
ma_steady_state_full_length = [oo_.dr.ys;oo_.dr.ys;oo_.dr.ys];
ma_steady_state             = ma_steady_state_full_length(MA_.iyv>0);

% evaluate derivatives at steady state
if options_.order==1
[waste,MA_.df1] = feval( [M_.fname '_dynamic'], ma_steady_state, [zeros(it_+1,M_.exo_nbr) ],  M_.params, oo_.dr.ys, it_);
elseif options_.order==2
[waste,MA_.df1,MA_.df2] = feval( [M_.fname '_dynamic'], ma_steady_state, [zeros(it_+1,M_.exo_nbr) ],  M_.params, oo_.dr.ys, it_);
elseif options_.order==3
[waste,MA_.df1,MA_.df2,MA_.df3] = feval( [M_.fname '_dynamic'], ma_steady_state, [zeros(it_+1,M_.exo_nbr) ],  M_.params, oo_.dr.ys, it_);
end
%sort the derivatives into dr order
temp=M_.lead_lag_incidence(:,oo_.dr.order_var);
temp=[temp(1,oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred) temp(2,:) temp(3,end-(oo_.dr.nfwrd+oo_.dr.nboth)+1:end)];
if M_.exo_nbr>0
    temp=[temp max(temp)+1:max(temp)+M_.exo_nbr];
end
temp=transpose(temp);
MA_.df1=MA_.df1(:,temp);
if options_.order>=2
J=temp(:,ones(1,length(temp))).';J=J(:);
I=(1:length(temp)).';I=I(:,ones(1,length(temp)));I=temp(I,1);
IJ=I+(J-1).*length(temp);
MA_.df2=MA_.df2(:,IJ);
if options_.order==3
K=temp(:,ones(1,length(temp)^2)).';K=K(:);
J=temp(:,ones(1,length(temp))).';J=J(:);I=(1:length(temp)^2).';I=I(:,ones(1,length(temp)));J=J(I,1);
I=(1:length(temp)).';I=I(:,ones(1,length(temp)^2));I=temp(I,1);
IJK=I+(J-1).*length(temp)+(K-1).*length(temp)^2;
MA_.df3=MA_.df3(:,IJK);
end
end
MA_.temp=temp;
%Start computing
disp('Computing First-Order Terms')
MA_=first_order(M_,MA_,oo_,options_);
disp('Done!')
if options_.order>=2
    disp('Computing Second-Order Terms')
    MA_=second_order(M_,MA_,oo_,options_);
    disp('Done!')
    if options_.order==3
        disp('Computing Third-Order Terms')
             MA_=third_order(M_,MA_,oo_,options_);
        disp('Done!')
     end
end
% 
if MA_.calculate_kernels==1
    disp('Computing Kernels. This may take a while.')
    [MA_]=kernels(M_,MA_,oo_,options_);
    disp('Done!')
    if MA_.plot_kernels==1
            MA_=plots_kernels(M_,MA_,oo_,options_);
     end
      
end
 if MA_.calculate_simulations==1
     disp('Computing Simulations.')
     if nargin==5
         [MA_]=MA_simul(M_,MA_,options_,oo_,varargin{2});
     else
         [MA_]=MA_simul(M_,MA_,options_,oo_,M_.endo_names);
     end
     disp('Done!')
 end
 if MA_.calculate_irf==1
         disp('Computing Impulse Responses')
         [MA_]=calculate_irfs(M_,MA_,options_,oo_);
          disp('Done!')
          if MA_.plot_irf==1
              if nargin==5
         [MA_]=plots_irf(M_,MA_,options_,oo_,varargin{2});
              else
         [MA_]=plots_irf(M_,MA_,options_,oo_,M_.endo_names);
             end
          end
 end
time_to_run=toc;
hours_to_run=floor(time_to_run/3600);
minutes_to_run=floor((time_to_run-hours_to_run*3600)/60);
seconds_to_run=time_to_run-hours_to_run*3600-minutes_to_run*60;
disp(['nonlinear_MA Execution Time:' blanks(2) num2str(hours_to_run) ' hour(s), ' num2str(minutes_to_run) ' minute(s), and '  num2str(seconds_to_run) ' second(s)'])
    