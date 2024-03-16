function [MA_]=calculate_irfs(M_,MA_,options_,oo_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate_irfs.m
%
% This file calulates IRF's recursively using the linearly recursive
% state-space representation associated with the nonlinear moving average
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

select_state=oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred;
MA_.irf.first=zeros(M_.exo_nbr*M_.endo_nbr,options_.irf);
MA_.irf.total=zeros(M_.exo_nbr*M_.endo_nbr,options_.irf);
if options_.order>=2
    MA_.irf.second=zeros(M_.exo_nbr*M_.endo_nbr,options_.irf);
if options_.order==3
    MA_.irf.third=zeros(M_.exo_nbr*M_.endo_nbr,options_.irf);
    MA_.irf.first_sigma_2=zeros(M_.exo_nbr*M_.endo_nbr,options_.irf);
end
end
for jj=1:M_.exo_nbr
S=zeros(M_.exo_nbr,1);S(jj,1)=1;
simul_first(:,1)=MA_.BETA_0*S;
if options_.order>=2
    exe=alt_kron(S,S);
simul_second(:,1)=(1/2)*MA_.BETA_00*exe;
if options_.order==3
 simul_first_sigma_2(:,1)=(1/2)*MA_.BETA_sigma_2_0*S;
simul_third(:,1)=(1/6)*MA_.BETA_000*alt_kron(S,exe);
end
end
for t=2:options_.irf
    simul_first(:,t)=MA_.ALPHA*simul_first(select_state,t-1);
if options_.order>=2
    SxS=alt_kron(simul_first(select_state,t-1),simul_first(select_state,t-1));
simul_second(:,t)=MA_.ALPHA*simul_second(select_state,t-1)+(1/2)*MA_.BETA_22*SxS;
if options_.order==3
    simul_first_sigma_2(:,t)=MA_.ALPHA*simul_first_sigma_2(select_state,t-1)+(1/2)*MA_.BETA_sigma_2_1*simul_first(select_state,t-1);
simul_third(:,t)=MA_.ALPHA*simul_third(select_state,t-1)+(1/6)*MA_.BETA_333_1*alt_kron(simul_first(select_state,t-1),SxS)...
    +MA_.BETA_22*alt_kron(simul_second(select_state,t-1),simul_first(select_state,t-1));
end
end
end
MA_.irf.first((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)=simul_first(oo_.dr.inv_order_var,:);
if options_.order>=2
    MA_.irf.second((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)=simul_second(oo_.dr.inv_order_var,:);
if options_.order==3
    MA_.irf.third((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)=simul_third(oo_.dr.inv_order_var,:);
    MA_.irf.first_sigma_2((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)=simul_first_sigma_2(oo_.dr.inv_order_var,:);
end
end
if options_.order==1
     MA_.irf.total((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)=MA_.irf.first((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2);
elseif options_.order==2
     MA_.irf.total((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)=MA_.irf.first((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2)+MA_.irf.second((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)*((MA_.shock_scale^2)*M_.Sigma_e(jj,jj));
elseif options_.order==3
    MA_.irf.total((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)=(MA_.irf.first((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)+MA_.irf.first_sigma_2((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:))*MA_.shock_scale*(M_.Sigma_e(jj,jj))^(1/2)+MA_.irf.second((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)*((MA_.shock_scale^2)*M_.Sigma_e(jj,jj))+MA_.irf.third((jj-1)*M_.endo_nbr+1:jj*M_.endo_nbr,:)*(MA_.shock_scale^3)*(M_.Sigma_e(jj,jj))^(3/2);
end

end

