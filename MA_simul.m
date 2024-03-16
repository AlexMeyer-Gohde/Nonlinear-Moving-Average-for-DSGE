 function [MA_]=MA_simul(M_,MA_,options_,oo_,var_list_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MA_simul.m
%
% This file calculates simulation(s) from a sequence of random shocks
% using the linearly recursive state-space representation associated with 
% the nonlinear moving average and compares the simulation with the one 
% produced by dynare in a plot
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

select_state=oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred;
%Endogenous pruning algorithm begins here
MA_.simul_first=zeros(M_.endo_nbr,options_.periods);
if options_.order>=2
    MA_.simul_second=zeros(M_.endo_nbr,options_.periods);
if options_.order==3
    MA_.simul_third=zeros(M_.endo_nbr,options_.periods);
    MA_.simul_first_sigma_2=zeros(M_.endo_nbr,options_.periods);
end
end
E=oo_.exo_simul';
MA_.simul_first(:,1)=MA_.BETA_0*E(:,1);
if options_.order>=2
    exe=alt_kron(E(:,1),E(:,1));
MA_.simul_second(:,1)=(1/2)*MA_.BETA_00*exe;
if options_.order==3
 MA_.simul_first_sigma_2(:,1)=(1/2)*MA_.BETA_sigma_2_0*E(:,1);
MA_.simul_third(:,1)=(1/6)*MA_.BETA_000*alt_kron(E(:,1),exe);

end
end
for t=2:options_.periods
    MA_.simul_first(:,t)=MA_.ALPHA*MA_.simul_first(select_state,t-1)+MA_.BETA_0*E(:,t);
if options_.order>=2
    exe=alt_kron(E(:,t),E(:,t));sxe=alt_kron(MA_.simul_first(select_state,t-1),E(:,t));sxs=alt_kron(MA_.simul_first(select_state,t-1),MA_.simul_first(select_state,t-1));
MA_.simul_second(:,t)=MA_.ALPHA*MA_.simul_second(select_state,t-1)...
    +(1/2)*(MA_.BETA_22*sxs+2*MA_.BETA_20*sxe+MA_.BETA_00*exe);
if options_.order==3
    MA_.simul_first_sigma_2(:,t)=MA_.ALPHA*MA_.simul_first_sigma_2(select_state,t-1)...
        +(1/2)*(MA_.BETA_sigma_2_0*E(:,t)+MA_.BETA_sigma_2_1*MA_.simul_first(select_state,t-1));
MA_.simul_third(:,t)=MA_.ALPHA*MA_.simul_third(select_state,t-1)...
    +(1/6)*(MA_.BETA_333_1*alt_kron(MA_.simul_first(select_state,t-1),sxs)+MA_.BETA_000*alt_kron(E(:,t),exe))...
    +(1/2)*(MA_.BETA_330_1*alt_kron(sxs,E(:,t))+MA_.BETA_300*alt_kron(MA_.simul_first(select_state,t-1),exe))...
    +MA_.BETA_22*alt_kron(MA_.simul_second(select_state,t-1),MA_.simul_first(select_state,t-1))+MA_.BETA_20*alt_kron(MA_.simul_second(select_state,t-1),E(:,t));
% MA_.simul_third(:,t)=MA_.ALPHA*MA_.simul_third(select_state,t-1)...
%     +(1/6)*(own_A_times_B_kronecker_C(MA_.BETA_333_1,MA_.simul_first(select_state,t-1),sxs,options_,MA_)+own_A_times_B_kronecker_C(MA_.BETA_000,E(:,t),exe,options_,MA_))...
%     +(1/2)*(own_A_times_B_kronecker_C(MA_.BETA_330_1,sxs,E(:,t),options_,MA_)+own_A_times_B_kronecker_C(MA_.BETA_300,MA_.simul_first(select_state,t-1),exe,options_,MA_))...
%     +MA_.BETA_22*alt_kron(MA_.simul_second(select_state,t-1),MA_.simul_first(select_state,t-1))+MA_.BETA_20*alt_kron(MA_.simul_second(select_state,t-1),E(:,t));
end
end
end
%Endogenous pruning algorithm ends here
if options_.order==1
    MA_.simul_first=MA_.simul_first(oo_.dr.inv_order_var,:);
     MA_.simul_total=MA_.simul_first+repmat(oo_.dr.ys,[1 options_.periods]);
elseif options_.order==2
    MA_.simul_first=MA_.simul_first(oo_.dr.inv_order_var,:);
    MA_.simul_second=MA_.simul_second(oo_.dr.inv_order_var,:);
    Y_sigma_2=MA_.Y_sigma_2(oo_.dr.inv_order_var,:);
    MA_.simul_total=MA_.simul_second+MA_.simul_first+repmat(oo_.dr.ys+0.5*Y_sigma_2,[1 options_.periods]);
elseif options_.order==3
    MA_.simul_first=MA_.simul_first(oo_.dr.inv_order_var,:);
    MA_.simul_second=MA_.simul_second(oo_.dr.inv_order_var,:);
    Y_sigma_2=MA_.Y_sigma_2(oo_.dr.inv_order_var,:);
    MA_.simul_first_sigma_2=MA_.simul_first_sigma_2(oo_.dr.inv_order_var,:);
    MA_.simul_third=MA_.simul_third(oo_.dr.inv_order_var,:);
MA_.simul_total=MA_.simul_third+MA_.simul_first_sigma_2+MA_.simul_second+MA_.simul_first+repmat(oo_.dr.ys+0.5*Y_sigma_2,[1 options_.periods]);
end
%Begin combining and plotting


for ii=1:length(variable_select)
    if options_.order==1
if MA_.plot_simulations==1
  figure('Units','characters','Position',[1 1 240 70]);
    plot(1+options_.drop:options_.periods,MA_.simul_total(variable_select(ii),1+options_.drop:end),'b')
    eval(sprintf('title(''Simulation of %s'')',M_.endo_names(variable_select(ii),:)))
     ylabel('Value');
    xlabel('Periods');
    legend({'First-Order Accurate'},'Location','Best')
end
    elseif options_.order==2
        if MA_.plot_simulations==1
        figure('Units','characters','Position',[1 1 240 70]);
    subplot(3,3,[1:6]);plot(1+options_.drop:options_.periods,MA_.simul_total(variable_select(ii),1+options_.drop:end),'b',1+options_.drop:options_.periods,MA_.simul_first(variable_select(ii),1+options_.drop:end)+oo_.dr.ys(variable_select(ii)),'b-.',1+options_.drop:options_.periods,oo_.endo_simul(variable_select(ii),1+options_.drop:options_.periods),'k--')
    eval(sprintf('title(''Simulation of %s'')',M_.endo_names(variable_select(ii),:)))
     ylabel('Value');
    xlabel('Periods');
    if options_.pruning==1
    legend({'Second-Order Accurate','First-Order Accurate','Dynare: Pruned Second Order'},'Location','Best')
    else
    legend({'Second-Order Accurate','First-Order Accurate','Dynare: Second Order'},'Location','Best')
    end
       subplot(3,3,[7]); plot(1+options_.drop:options_.periods,MA_.simul_first(variable_select(ii),1+options_.drop:end),'b')
    title('First-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot(1+options_.drop:options_.periods,repmat(0, [1,options_.periods-options_.drop]),'k-')
     hold off
      subplot(3,3,[8]); plot(1+options_.drop:options_.periods,MA_.simul_second(variable_select(ii),1+options_.drop:end),'b')
    title('Second-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot(1+options_.drop:options_.periods,repmat(0, [1,options_.periods-options_.drop]),'k-')
     hold off
      subplot(3,3,[9]); plot(1+options_.drop:options_.periods,repmat(oo_.dr.ys(variable_select(ii)), [1 options_.periods-options_.drop]),'b',1+options_.drop:options_.periods,repmat(oo_.dr.ys(variable_select(ii))+0.5*Y_sigma_2(variable_select(ii)), [1 options_.periods-options_.drop]),'b-.')
    title('Steady-State and Constant')
    legend({'Steady-State','plus Risk Adjutsment'},'Location','Best')  
        end
    elseif options_.order==3
                if MA_.plot_simulations==1
        figure('Units','characters','Position',[1 1 240 70]);
     subplot(4,3,[1:6]);plot(1+options_.drop:options_.periods,MA_.simul_total(variable_select(ii),1+options_.drop:end),'b',1+options_.drop:options_.periods,oo_.dr.ys(variable_select(ii))+0.5*Y_sigma_2(variable_select(ii))+MA_.simul_second(variable_select(ii),1+options_.drop:end)+MA_.simul_first(variable_select(ii),1+options_.drop:end),'b--',1+options_.drop:options_.periods,MA_.simul_first(variable_select(ii),1+options_.drop:end)+oo_.dr.ys(variable_select(ii)),'b-.',1+options_.drop:options_.periods,oo_.endo_simul(variable_select(ii),1+options_.drop:options_.periods),'k--')
    eval(sprintf('title(''Simulation of %s'')',M_.endo_names(variable_select(ii),:)))
     ylabel('Value');
    xlabel('Periods');
    legend({'Third-Order Accurate','Second-Order Accurate','First-Order Accurate','Dynare: Third Order'},'Location','Best')
       subplot(4,3,[7]); plot(1+options_.drop:options_.periods,MA_.simul_first(variable_select(ii),1+options_.drop:end),'b')
    title('First-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([1+options_.drop:options_.periods],repmat(0, [1,options_.periods-options_.drop]),'k-')
     hold off
      subplot(4,3,[8]); plot(1+options_.drop:options_.periods,MA_.simul_second(variable_select(ii),1+options_.drop:end),'b')
    title('Second-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([1+options_.drop:options_.periods],repmat(0, [1,options_.periods-options_.drop]),'k-')
     hold off
     subplot(4,3,[10]); plot(1+options_.drop:options_.periods,MA_.simul_first_sigma_2(variable_select(ii),1+options_.drop:end),'b')
    title('Risk Correction to First-Order')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([1+options_.drop:options_.periods],repmat(0, [1,options_.periods-options_.drop]),'k-')
     hold off
      subplot(4,3,[11]); plot(1+options_.drop:options_.periods,MA_.simul_third(variable_select(ii),1+options_.drop:end),'b')
    title('Third-Order Component')
    ylabel('Deviations');
    xlabel('Periods');
    hold on
     plot([1+options_.drop:options_.periods],repmat(0, [1,options_.periods-options_.drop]),'k-')
     hold off
      subplot(4,3,[9]); plot(1+options_.drop:options_.periods,repmat(oo_.dr.ys(variable_select(ii)), [1 options_.periods-options_.drop]),'b',1+options_.drop:options_.periods,repmat(oo_.dr.ys(variable_select(ii))+0.5*Y_sigma_2(variable_select(ii)), [1 options_.periods-options_.drop]),'b-.')
    title('Steady-State and Constant')
    legend({'Steady-State','plus Risk Adjutsment'},'Location','Best')   
                end
    end
end
