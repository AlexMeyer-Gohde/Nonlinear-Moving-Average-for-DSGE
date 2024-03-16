function [MA_]=third_order(M_,MA_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% third_order.m
%
% This file produces the coefficients for a third-order approximation in 
% addition to those produced by first_order.m and second_order.m 
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

%Beta_333_1
[temp]  = own_A_times_B_kronecker_C(MA_.BETA_22(MA_.BETA_0_select_fwdendo,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),options_,MA_);
MA_.GAMMA_22=[sparse(oo_.dr.npred,oo_.dr.npred^2);MA_.BETA_22;MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,:)*MA_.BETA_22(MA_.BETA_0_select_state,:)+temp;sparse(M_.exo_nbr,oo_.dr.npred^2)];
%syl_DD=sparse_kron_prod(MA_.df2,MA_.GAMMA_22,MA_.GAMMA_1);
[syl_DD]=own_sparse_hessian_times_B_kronecker_C(MA_.df2,full(MA_.GAMMA_22),MA_.GAMMA_1,options_,MA_);
[temp]  = own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22(MA_.BETA_0_select_fwdendo,:),MA_.BETA_22(MA_.BETA_0_select_state,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),options_,MA_);
syl_DD=(syl_DD+temp)*(speye(oo_.dr.npred^3)+alt_kron(speye(oo_.dr.npred),commutation_sparse(oo_.dr.npred,oo_.dr.npred))+commutation_sparse(oo_.dr.npred^2,oo_.dr.npred));
temp=sparse_kron_prod_3(MA_.df3,MA_.GAMMA_1,MA_.GAMMA_1,MA_.GAMMA_1);
syl_DD=syl_DD+temp;
[waste, MA_.BETA_333_1_FWDENDO] = gensylv(3,MA_.AA_0*MA_.syl_BB,MA_.AA_0*MA_.syl_CC,MA_.ALPHA(MA_.ALPHA_row_select_state,:),-MA_.AA_0*syl_DD);
[temp]  = own_A_times_B_kronecker_C(MA_.syl_CC*MA_.BETA_333_1_FWDENDO,MA_.ALPHA(MA_.ALPHA_row_select_state,:),alt_kron(MA_.ALPHA(MA_.ALPHA_row_select_state,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:)),options_,MA_);
MA_.BETA_333_1_OTHER=- MA_.AA_pinv*(MA_.syl_BB*MA_.BETA_333_1_FWDENDO+temp+syl_DD);
MA_.BETA_333_1=[MA_.BETA_333_1_OTHER;MA_.BETA_333_1_FWDENDO];

%Beta_330_1
[temp]  = own_A_times_B_kronecker_C(MA_.BETA_22(MA_.BETA_0_select_fwdendo,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:),options_,MA_);
MA_.GAMMA_20=[sparse(oo_.dr.npred,M_.exo_nbr*oo_.dr.npred);MA_.BETA_20;MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,:)*MA_.BETA_20(MA_.BETA_0_select_state,:)+temp;sparse(M_.exo_nbr,M_.exo_nbr*oo_.dr.npred)];
[syl_DD ]  = own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_333_1(MA_.BETA_0_select_fwdendo,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),alt_kron(MA_.ALPHA(MA_.ALPHA_row_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:)),options_,MA_);
[temp]  = own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22(MA_.BETA_0_select_fwdendo,:),MA_.BETA_20(MA_.ALPHA_row_select_state,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),options_,MA_);
syl_DD=syl_DD+temp*(alt_kron(speye(oo_.dr.npred),commutation_sparse(M_.exo_nbr,oo_.dr.npred))+commutation_sparse(oo_.dr.npred*M_.exo_nbr,oo_.dr.npred));
[temp]  = own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22(MA_.BETA_0_select_fwdendo,:),MA_.BETA_22(MA_.ALPHA_row_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:),options_,MA_);
syl_DD=syl_DD+temp;
temp=sparse_kron_prod_3(MA_.df3,MA_.GAMMA_1,MA_.GAMMA_1,MA_.X_0);
syl_DD=syl_DD+temp;
%temp=sparse_kron_prod(MA_.df2,MA_.GAMMA_1,MA_.GAMMA_20);
[temp]  = own_sparse_hessian_times_B_kronecker_C(MA_.df2,MA_.GAMMA_1,full(MA_.GAMMA_20),options_,MA_);
%syl_DD=syl_DD+temp*(speye(M_.exo_nbr*oo_.dr.npred^2)+alt_kron(speye(oo_.dr.npred),commutation_sparse(M_.exo_nbr,oo_.dr.npred))); HL->AMG:Replace this line with the line below, 28 March 2013
syl_DD=syl_DD+temp*(speye(M_.exo_nbr*oo_.dr.npred^2)+alt_kron(commutation_sparse(oo_.dr.npred,oo_.dr.npred),speye(M_.exo_nbr)) );
%temp=sparse_kron_prod(MA_.df2,MA_.GAMMA_22,MA_.X_0);
[temp]  = own_sparse_hessian_times_B_kronecker_C(MA_.df2,full(MA_.GAMMA_22),MA_.X_0,options_,MA_);
syl_DD=syl_DD+temp;
MA_.BETA_330_1=-MA_.syl_AA_tilde\syl_DD;

%Beta_300
[syl_DD]= own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22_FWDENDO,MA_.BETA_20(MA_.BETA_0_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:),options_,MA_);
syl_DD=syl_DD*(speye(oo_.dr.npred*M_.exo_nbr^2)+alt_kron(speye(oo_.dr.npred),commutation_sparse(M_.exo_nbr,M_.exo_nbr)));
[temp]= own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22_FWDENDO,MA_.BETA_00(MA_.BETA_0_select_state,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),options_,MA_);
syl_DD=syl_DD+temp*commutation_sparse(M_.exo_nbr^2,oo_.dr.npred);
[temp]  = own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_333_1(MA_.BETA_0_select_fwdendo,:),MA_.ALPHA(MA_.ALPHA_row_select_state,:),alt_kron(MA_.BETA_0(MA_.BETA_0_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:)),options_,MA_);
syl_DD=syl_DD+temp;
%temp= sparse_kron_prod(MA_.df2,MA_.GAMMA_20,MA_.X_0);
[temp]  = own_sparse_hessian_times_B_kronecker_C(MA_.df2,full(MA_.GAMMA_20),MA_.X_0,options_,MA_);
syl_DD=syl_DD+temp*(speye(oo_.dr.npred*M_.exo_nbr^2)+alt_kron(speye(oo_.dr.npred),commutation_sparse(M_.exo_nbr,M_.exo_nbr)));
MA_.X_00 = [ sparse(oo_.dr.npred,M_.exo_nbr^2);MA_.BETA_00; MA_.ALPHA(MA_.ALPHA_row_select_fwdendo,:)*MA_.BETA_00(MA_.BETA_0_select_state,:)+MA_.BETA_22(MA_.BETA_0_select_fwdendo,:)*alt_kron(MA_.BETA_0(MA_.BETA_0_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:));sparse(M_.exo_nbr,M_.exo_nbr^2)];
%temp= sparse_kron_prod(MA_.df2,MA_.GAMMA_1,MA_.X_00);
[temp]  = own_sparse_hessian_times_B_kronecker_C(MA_.df2,MA_.GAMMA_1,full(MA_.X_00),options_,MA_);
syl_DD=syl_DD+temp;
temp= sparse_kron_prod_3(MA_.df3,MA_.GAMMA_1,MA_.X_0,MA_.X_0);
syl_DD=syl_DD+temp;
MA_.BETA_300=-MA_.syl_AA_tilde\syl_DD;

%Beta_000
[syl_DD]= own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_22_FWDENDO,MA_.BETA_00(MA_.BETA_0_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:),options_,MA_);
syl_DD=syl_DD*(speye(M_.exo_nbr^3)+alt_kron(speye(M_.exo_nbr),commutation_sparse(M_.exo_nbr,M_.exo_nbr))+commutation_sparse(M_.exo_nbr^2,M_.exo_nbr));
[temp]  = own_A_times_B_kronecker_C(MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_333_1(MA_.BETA_0_select_fwdendo,:),MA_.BETA_0(MA_.BETA_0_select_state,:),alt_kron(MA_.BETA_0(MA_.BETA_0_select_state,:),MA_.BETA_0(MA_.BETA_0_select_state,:)),options_,MA_);
syl_DD=syl_DD+temp;
%temp= sparse_kron_prod(MA_.df2,MA_.X_00,MA_.X_0);
[temp]  = own_sparse_hessian_times_B_kronecker_C(MA_.df2,full(MA_.X_00),MA_.X_0,options_,MA_);
syl_DD=syl_DD+temp*(speye(M_.exo_nbr^3)+alt_kron(speye(M_.exo_nbr),commutation_sparse(M_.exo_nbr,M_.exo_nbr))+commutation_sparse(M_.exo_nbr^2,M_.exo_nbr));
temp= sparse_kron_prod_3(MA_.df3,MA_.X_0,MA_.X_0,MA_.X_0);
syl_DD=syl_DD+temp;
MA_.BETA_000=-MA_.syl_AA_tilde\syl_DD;
 
%Now the time-varying uncertainty correction terms
 %Selection indices used to extract relevant components of tensor
 %derivatives
 temp_1=MA_.df1_select_fwdendo_tplus;temp_1=temp_1';l_1=length(temp_1);
 temp_2=1:oo_.dr.npred+M_.endo_nbr+oo_.dr.nboth+oo_.dr.nfwrd+M_.exo_nbr;temp_2=temp_2';l_2=length(temp_2);
J=temp_1(:,ones(1,l_2)).';J=J(:);
I=(1:l_2).';I=I(:,ones(1,l_1));I=temp_2(I,1);
MA_.df2_select_fwdendo_tplus_x=I+(J-1).*l_2;

K=temp_1(:,ones(1,l_2*l_1)).';K=K(:);
J=temp_1(:,ones(1,l_2)).';J=J(:);JJ=(1:l_1*l_2).';JJ=JJ(:,ones(1,l_1));J=J(JJ,1);
I=(1:l_2).';I=I(:,ones(1,l_1^2));I=temp_2(I,1);
MA_.df3_select_fwdendo_tplus_fwdendo_tplus_x=I+(K-1).*(l_2*l_2)+(J-1).*l_2;

%Beta_sigma_2_1
syl_DD =sparse_kron_prod_3(MA_.df3(:,MA_.df3_select_fwdendo_tplus_fwdendo_tplus_x),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.GAMMA_1);
%temp=sparse_kron_prod(MA_.df2(:,MA_.df2_select_fwdendo_tplus_x),MA_.BETA_00(MA_.BETA_0_select_fwdendo,:),MA_.GAMMA_1);
[temp]  = own_A_times_B_kronecker_C(full(MA_.df2(:,MA_.df2_select_fwdendo_tplus_x)),MA_.BETA_00(MA_.BETA_0_select_fwdendo,:),MA_.GAMMA_1,options_,MA_);
syl_DD=syl_DD+temp;
%temp  = sparse_kron_prod(2*MA_.df2(:,MA_.df2_select_fwdendo_tplus_fwdendo_tplus),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.BETA_20(MA_.BETA_0_select_fwdendo,:)*commutation_sparse(oo_.dr.npred,M_.exo_nbr));
[temp]  = own_sparse_hessian_times_B_kronecker_C(2*(MA_.df2(:,MA_.df2_select_fwdendo_tplus_fwdendo_tplus)),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.BETA_20(MA_.BETA_0_select_fwdendo,:)*commutation(oo_.dr.npred,M_.exo_nbr),options_,MA_);%%%AMG 23.01.13 removed commutation%%%%%
temp= temp+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_300(MA_.BETA_0_select_fwdendo,:)*commutation_sparse(oo_.dr.npred,M_.exo_nbr^2);
syl_DD=syl_DD+temp*alt_kron(speye(M_.exo_nbr^2),MA_.ALPHA(MA_.ALPHA_row_select_state,:));
%syl_DD=sparse_kron_prod(syl_DD,M_.Sigma_e(:),speye(oo_.dr.npred));
[syl_DD]  = own_A_times_B_kronecker_C(full(syl_DD),M_.Sigma_e(:),eye(oo_.dr.npred),options_,MA_);
%temp  = sparse_kron_prod(MA_.df2,[MA_.Y_sigma_2(MA_.BETA_0_select_state);MA_.Y_sigma_2;MA_.Y_sigma_2(MA_.BETA_0_select_fwdendo);zeros(M_.exo_nbr,1)],MA_.GAMMA_1);
[temp]  = own_sparse_hessian_times_B_kronecker_C(MA_.df2,[MA_.Y_sigma_2(MA_.BETA_0_select_state);MA_.Y_sigma_2;MA_.Y_sigma_2(MA_.BETA_0_select_fwdendo);zeros(M_.exo_nbr,1)],MA_.GAMMA_1,options_,MA_);
syl_DD=syl_DD+temp;
[waste, MA_.BETA_sigma_2_1_FWDENDO] = gensylv(1,MA_.AA_0*MA_.syl_BB,MA_.AA_0*MA_.syl_CC,MA_.ALPHA(MA_.ALPHA_row_select_state,:),-MA_.AA_0*syl_DD);
MA_.BETA_sigma_2_1_OTHER=- MA_.AA_pinv*(MA_.syl_BB*MA_.BETA_sigma_2_1_FWDENDO+MA_.syl_CC*MA_.BETA_sigma_2_1_FWDENDO*MA_.ALPHA(MA_.ALPHA_row_select_state,:)+syl_DD);
MA_.BETA_sigma_2_1=[MA_.BETA_sigma_2_1_OTHER;MA_.BETA_sigma_2_1_FWDENDO];

%Beta_sigma_2_0
syl_DD =sparse_kron_prod_3(MA_.df3(:,MA_.df3_select_fwdendo_tplus_fwdendo_tplus_x),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.X_0);
%temp=sparse_kron_prod(MA_.df2(:,MA_.df2_select_fwdendo_tplus_x),MA_.BETA_00(MA_.BETA_0_select_fwdendo,:),MA_.X_0);
[temp]  = own_A_times_B_kronecker_C(full(MA_.df2(:,MA_.df2_select_fwdendo_tplus_x)),MA_.BETA_00(MA_.BETA_0_select_fwdendo,:),MA_.X_0,options_,MA_);
syl_DD=syl_DD+temp;
%temp  = sparse_kron_prod(2*MA_.df2(:,MA_.df2_select_fwdendo_tplus_fwdendo_tplus),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.BETA_20(MA_.BETA_0_select_fwdendo,:)*commutation_sparse(oo_.dr.npred,M_.exo_nbr)*alt_kron(speye(M_.exo_nbr),MA_.BETA_0(MA_.BETA_0_select_state,:)));
[temp]  = own_sparse_hessian_times_B_kronecker_C(2*(MA_.df2(:,MA_.df2_select_fwdendo_tplus_fwdendo_tplus)),MA_.BETA_0(MA_.BETA_0_select_fwdendo,:),MA_.BETA_20(MA_.BETA_0_select_fwdendo,:)*commutation(oo_.dr.npred,M_.exo_nbr)*alt_kron(eye(M_.exo_nbr),MA_.BETA_0(MA_.BETA_0_select_state,:)),options_,MA_);
syl_DD=syl_DD+temp;
temp= MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_300(MA_.BETA_0_select_fwdendo,:)*commutation_sparse(oo_.dr.npred,M_.exo_nbr^2)*alt_kron(speye(M_.exo_nbr^2),MA_.BETA_0(MA_.BETA_0_select_state,:));
syl_DD=syl_DD+temp;
[syl_DD]  = own_A_times_B_kronecker_C(syl_DD,M_.Sigma_e(:),eye(M_.exo_nbr),options_,MA_);
syl_DD  = syl_DD+MA_.df1(:,MA_.df1_select_fwdendo_tplus)*MA_.BETA_sigma_2_1_FWDENDO*MA_.BETA_0(MA_.BETA_0_select_state,:);
%temp  = sparse_kron_prod(MA_.df2,[MA_.Y_sigma_2(MA_.BETA_0_select_state);MA_.Y_sigma_2;MA_.Y_sigma_2(MA_.BETA_0_select_fwdendo);zeros(M_.exo_nbr,1)],MA_.X_0);
[temp]  = own_sparse_hessian_times_B_kronecker_C(MA_.df2,[MA_.Y_sigma_2(MA_.BETA_0_select_state);MA_.Y_sigma_2;MA_.Y_sigma_2(MA_.BETA_0_select_fwdendo);zeros(M_.exo_nbr,1)],MA_.X_0,options_,MA_);
syl_DD=syl_DD+temp;
MA_.BETA_sigma_2_0=-MA_.syl_AA_tilde\syl_DD;
