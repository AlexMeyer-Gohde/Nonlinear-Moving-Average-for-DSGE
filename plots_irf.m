 function [MA_]=plots_irf(M_,MA_,options_,oo_,var_list_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_irf.m
%
% This file produces plots of the IRF's calculated in calculate_irfs,
% including a decomposition of the IRF into contributing orders of
% nonlinearity and uncertainty.
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


if isempty(var_list_)==0
[waste, variable_select] = ismember(cellstr(var_list_), cellstr(M_.endo_names));
else
    variable_select=1:M_.endo_nbr;
end
if options_.order>=2
    Y_sigma_2=MA_.Y_sigma_2(oo_.dr.inv_order_var,:);
end

irf_length_vector=1:options_.irf;
for jj=1:M_.exo_nbr
for ii=1:length(variable_select)
    
    if options_.order==1
  figure('Units','characters','Position',[1 1 240 70]);
    plot(irf_length_vector-1,MA_.irf.total((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector),'b')
    eval(sprintf('title(''Impulse Response of %s to a %s Std. Dev. Shock in %s'')',M_.endo_names(variable_select(ii),:),num2str((MA_.shock_scale)),M_.exo_names(jj,:)))
     ylabel('Deviations');
    xlabel('Periods since Shock Realization');
    legend({'First-Order Accurate'})
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
    
    elseif options_.order==2
    figure('Units','characters','Position',[1 1 240 70]);
    subplot(3,3,[1:6]); plot([irf_length_vector-1],MA_.irf.total((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector),'b',[irf_length_vector-1],MA_.irf.first((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2),'b-.')
    eval(sprintf('title(''Impulse Response of %s to a %s Std. Dev. Shock in %s'')',M_.endo_names(variable_select(ii),:),num2str((MA_.shock_scale)),M_.exo_names(jj,:)))
    ylabel('Deviations');
    xlabel('Periods since Shock Realization');
    legend({'Second-Order Accurate','First-Order Accurate'})
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
     subplot(3,3,[7]); plot([irf_length_vector-1],MA_.irf.first((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2),'b')
    title('First-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
      subplot(3,3,[8]); plot([irf_length_vector-1],MA_.irf.second((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*((MA_.shock_scale^2)*M_.Sigma_e(jj,jj)),'b')
    title('Second-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
      subplot(3,3,[9]); plot([irf_length_vector-1],repmat(oo_.dr.ys(variable_select(ii)), [1 options_.irf]),'b',[irf_length_vector-1],repmat(oo_.dr.ys(variable_select(ii))+Y_sigma_2(variable_select(ii)), [1 options_.irf]),'b-.')
    title('Steady-State and Constant')
    legend({'Steady-State','plus Risk Adjutsment'})
    
    elseif options_.order==3
    figure('Units','characters','Position',[1 1 240 70]);
    subplot(4,3,[1:6]); plot([irf_length_vector-1],MA_.irf.total((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector),'b',[irf_length_vector-1],MA_.irf.first((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2)+MA_.irf.second((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*((MA_.shock_scale^2)*M_.Sigma_e(jj,jj)),'b--',[irf_length_vector-1],MA_.irf.first((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2),'b-.')
    eval(sprintf('title(''Impulse Response of %s to a %s Std. Dev. Shock in %s'')',M_.endo_names(variable_select(ii),:),num2str((MA_.shock_scale)),M_.exo_names(jj,:)))
    ylabel('Deviations');
    xlabel('Periods since Shock Realization');
    legend({'Third-Order Accurate','Second-Order Accurate','First-Order Accurate'})
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
     subplot(4,3,[7]); plot([irf_length_vector-1],MA_.irf.first((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2),'b')
    title('First-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
      subplot(4,3,[8]); plot([irf_length_vector-1],MA_.irf.second((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*((MA_.shock_scale^2)*M_.Sigma_e(jj,jj)),'b')
    title('Second-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
          subplot(4,3,[10]); plot([irf_length_vector-1],MA_.irf.first_sigma_2((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2),'b')
    title('Risk Correction to First-Order')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
    set(gca,'xlim',[0 options_.irf])
      subplot(4,3,[11]); plot([irf_length_vector-1],MA_.irf.third((jj-1)*M_.endo_nbr+variable_select(ii),irf_length_vector)*(MA_.shock_scale^3)*(M_.Sigma_e(jj,jj))^(3/2),'b')
    title('Third-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([irf_length_vector-1],zeros(1, options_.irf),'k-')
     hold off
     set(gca,'xlim',[0 options_.irf])
      subplot(4,3,[9]); plot([irf_length_vector-1],repmat(oo_.dr.ys(variable_select(ii)), [1 options_.irf]),'b',[irf_length_vector-1],repmat(oo_.dr.ys(variable_select(ii))+Y_sigma_2(variable_select(ii)), [1 options_.irf]),'b-.')
    title('Steady-State and Constant')
    legend({'Steady-State','plus Risk Adjustment'})
    set(gca,'xlim',[0 options_.irf])
    end
end
end