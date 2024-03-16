function [MA_]=second_order(M_,MA_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% second_order.m
%
% This file produces the coefficients for a second-order approximation in 
% addition to those produced by first_order.m
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
   % build up lhs coefficients for Sylvesters
   MA_.df1_select_both_tminus=1+oo_.dr.npred-oo_.dr.nboth:oo_.dr.npred;
   MA_.df1_select_both_t=1+2*oo_.dr.npred-oo_.dr.nboth+oo_.dr.nstatic:2*oo_.dr.npred+oo_.dr.nstatic;
   MA_.df1_select_both_tplus=1+oo_.dr.npred+M_.endo_nbr:oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth;
   MA_.df1_select_fwdendo_tplus=1+oo_.dr.npred+M_.endo_nbr:oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd;
   MA_.ALPHA_row_select_fwdendo=1+oo_.dr.nstatic+oo_.dr.npred-oo_.dr.nboth:M_.endo_nbr;
   MA_.ALPHA_col_select_both=oo_.dr.npred-oo_.dr.nboth+1:oo_.dr.npred;
   MA_.df1_select_fwd_t=1+oo_.dr.npred+M_.endo_nbr-oo_.dr.nfwrd:oo_.dr.npred+M_.endo_nbr;
   MA_.df1_select_fwd_tplus=1+oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth:oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd;
   MA_.df1_select_st_t=1+oo_.dr.npred:oo_.dr.nstatic+oo_.dr.npred;
   MA_.df1_select_bwd_tminus=1:oo_.dr.npred-oo_.dr.nboth;
   MA_.df1_select_bwd_t=1+oo_.dr.npred+oo_.dr.nstatic:2*oo_.dr.npred+oo_.dr.nstatic-oo_.dr.nboth;
   MA_.df1_select_state_t= 1+oo_.dr.npred+oo_.dr.nstatic:2*oo_.dr.npred+oo_.dr.nstatic;
   MA_.ALPHA_col_select_bwd=1:oo_.dr.npred-oo_.dr.nboth;
   MA_.ALPHA_row_select_state=oo_.dr.nstatic+1:oo_.dr.npred+oo_.dr.nstatic;
   MA_.BETA_0_select_state=oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred;
   MA_.BETA_0_select_fwdendo=oo_.dr.nstatic+oo_.dr.npred-oo_.dr.nboth+1:M_.endo_nbr;
   
   MA_.syl_BB = [MA_.df1(:,MA_.df1_select_both_t)+ MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,MA_.ALPHA_col_select_both),...
                    MA_.df1(:,MA_.df1_select_fwd_t)];
     MA_.syl_CC = MA_.df1(:,MA_.df1_select_fwdendo_tplus);

        MA_.AA_0 = (null([MA_.df1(:,MA_.df1_select_st_t) MA_.df1(:,MA_.df1_select_bwd_t)+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,MA_.ALPHA_col_select_bwd)]'))';
        MA_.AA_pinv = pinv([MA_.df1(:,MA_.df1_select_st_t) MA_.df1(:,MA_.df1_select_bwd_t)+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,MA_.ALPHA_col_select_bwd)]);
        
% build up shifting matrix gamma(1), where (1) indicates that this shifting
% matrix bridges first order solution to the second order. 
% build up gamma(1), such that
%[y_(i-1)',y_(i)',y_(i+1)',U_(i)']'=gamma(1)*ys_(i-1). 
MA_.GAMMA_1 = [ eye(oo_.dr.npred)           
                                MA_.ALPHA                       
                                MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,:)*MA_.ALPHA(MA_.ALPHA_row_select_state,:);                    
                       zeros(M_.exo_nbr,oo_.dr.npred)             ];

% build up gamma(0), such that
%[y_(-1)',y_(0)',y_(1)',U_(0)']'
MA_.X_0 = [ zeros(oo_.dr.npred,M_.exo_nbr)           
                                MA_.BETA_0                     
                                MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,:)*MA_.BETA_0(MA_.BETA_0_select_state,:);                    
                       eye(M_.exo_nbr)             ];

        


% build the remaining second order-dependent matrices in the Sylvester
% temp=1:oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd;
% temp=temp';
% J=temp(:,ones(1,length(temp))).';J=J(:);
% I=(1:length(temp)).';I=I(:,ones(1,length(temp)));I=temp(I,1);
% IJ=I+(J-1).*(oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd+M_.exo_nbr);
% syl_DD =MA_.df2(:,IJ)*kron(MA_.GAMMA_1(1:oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd,:),MA_.GAMMA_1(1:oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd,:));
% syl_DD =MA_.df2*kron(MA_.GAMMA_1,MA_.GAMMA_1);
%syl_DD =sparse_kron_prod(MA_.df2,MA_.GAMMA_1,MA_.GAMMA_1);
[syl_DD] =own_sparse_hessian_times_B_kronecker_C(MA_.df2,MA_.GAMMA_1,MA_.GAMMA_1,options_,MA_);



% Call Dynare's Sylvester caller
[waste, MA_.BETA_22_FWDENDO] = gensylv(2,MA_.AA_0*MA_.syl_BB,MA_.AA_0*MA_.syl_CC,MA_.ALPHA(MA_.ALPHA_row_select_state,:),-MA_.AA_0*syl_DD);%can take advantage of Kronecker input on C in Dynare's Sylvester solver
%Zero out numbers numerically identical to zero according to the conditioning criterion

MA_.BETA_22_OTHER=- MA_.AA_pinv*(MA_.syl_BB*MA_.BETA_22_FWDENDO+MA_.syl_CC*MA_.BETA_22_FWDENDO*alt_kron(MA_.ALPHA(MA_.ALPHA_row_select_state,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:))+syl_DD);
MA_.BETA_22=[MA_.BETA_22_OTHER;MA_.BETA_22_FWDENDO];

MA_.syl_AA_tilde=[MA_.df1(:,MA_.df1_select_st_t),MA_.df1(:,MA_.df1_select_state_t)+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,:),MA_.df1(:,MA_.df1_select_fwd_t)];

[temp] =own_sparse_hessian_times_B_kronecker_C(MA_.df2,MA_.GAMMA_1,MA_.X_0,options_,MA_);
MA_.BETA_20=-MA_.syl_AA_tilde\(temp+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22_FWDENDO*alt_kron(MA_.ALPHA(MA_.ALPHA_row_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:)));

MA_.BETA_00=-MA_.syl_AA_tilde \(MA_.df2*alt_kron(MA_.X_0,MA_.X_0)+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22_FWDENDO*alt_kron(MA_.BETA_0(MA_.BETA_0_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:)));


% %This builds the term Y_{\tilde{epsilon}^2} using the second-order
% %solution. It measure the second-order dependence of the system's variables to a shock
% %in t+1. Only y_{t+1} has such a dependence.
% MA_.Y_eps_2=[zeros(2*M_.endo_nbr,M_.exo_nbr^2);MA_.BETA_2*kron([zeros(M_.endo_nbr,M_.exo_nbr);eye(M_.exo_nbr,M_.exo_nbr)],[zeros(M_.endo_nbr,M_.exo_nbr);eye(M_.exo_nbr,M_.exo_nbr)]);zeros(M_.exo_nbr,M_.exo_nbr^2)];  
% 
% %Now calculate Y_sigma_sigma
% MA_.Y_sigma_2=  -(2*MA_.matrix_polynomial_one)\((MA_.df2*kron(MA_.g_picker*MA_.Y_eps_1,MA_.g_picker*MA_.Y_eps_1)+MA_.df1*MA_.g_picker*MA_.Y_eps_2)*vec(M_.Sigma_e));
% %Zero out numbers numerically identical to zero according to the conditioning criterion
% MA_.Y_sigma_2(abs(MA_.Y_sigma_2)<MA_.condn)=0;
MA_.matrix_polynomial_one=[MA_.df1(:,MA_.df1_select_st_t), ...
    MA_.df1(:,MA_.df1_select_bwd_tminus)+MA_.df1(:,MA_.df1_select_bwd_t),...
    MA_.df1(:,MA_.df1_select_both_tminus)+MA_.df1(:,MA_.df1_select_both_t)+MA_.df1(:,MA_.df1_select_both_tplus),...
    MA_.df1(:,MA_.df1_select_fwd_t)+MA_.df1(:,MA_.df1_select_fwd_tplus)];


temp=MA_.df1_select_fwdendo_tplus;
temp=temp';
J=temp(:,ones(1,length(temp))).';J=J(:);
I=(1:length(temp)).';I=I(:,ones(1,length(temp)));I=temp(I,1);
MA_.df2_select_fwdendo_tplus_fwdendo_tplus=I+(J-1).*(oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd+M_.exo_nbr);


MA_.Y_sigma_2=  -(MA_.matrix_polynomial_one)\((MA_.df2(:,MA_.df2_select_fwdendo_tplus_fwdendo_tplus)*kron(MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:))+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_00(MA_.BETA_0_select_fwdendo,:))*M_.Sigma_e(:));