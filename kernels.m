 function [MA_]=kernels(M_,MA_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% kernels.m
%
% This file calculates the Volterra kernels
% 
%Included in the package nonlinear_MA
%
%THIS VERSION: 1.0.1 Dec. 19, 2011
%
%Copyright: Hong Lan and Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



select_state=oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred;
%Endogenous pruning algorithm begins here

MA_.max_length=MA_.max_kernel_length;
 MA_.Y_i=cell(1,MA_.max_length);
%Compute the first-order kernel
 MA_.Y_i{1}=MA_.BETA_0;
for i=2:MA_.max_length
    MA_.Y_i{i}=MA_.ALPHA*MA_.Y_i{i-1}(select_state,:);
end
%Compute the second-order kernel
if options_.order >= 2
%     %preallocate
    commute_ji_ij=commutation(M_.exo_nbr, M_.exo_nbr);
%     D=cell2mat(MA_.Y_i);
%     D=D(select_state,:);
%     D=mat2cell(D,[oo_.dr.npred],repmat(M_.exo_nbr,[1 MA_.max_length]));
%     D=reshape(D,[MA_.max_length,1]);%column vector of cells
%     D=[D, repmat({zeros(M_.exo_nbr,M_.exo_nbr)},[MA_.max_length 1] )]; %add a column of cells adjacent with each cell being a matrix of zeros (response of shock to past shock is zero with no serial correlation) 
%      D = reshape( D, [ 1 MA_.max_length 1 2 ] );%from here
%     D = permute( D, [ 1 4 3 2 ] );
%     D = reshape( D, [ 2 MA_.max_length] );%to here performs a block transpose
%     D=reshape(D,[2*MA_.max_length,1]);%here I vectorize columnwise
%     D=cell2mat(D);%and convert the column vector of cells into a matrix
%     DkronD=alt_kron(D,D);%kronecker the two together. This gives all possible combinations of: kron([MA_.Y_i{j-1};zeros(M_.exo_nbr,M_.exo_nbr)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)]);
%     select=vec(repmat([0:1:(oo_.dr.npred)-1]'*MA_.max_length*(oo_.dr.npred),[1,(oo_.dr.npred)])')+repmat([1:(oo_.dr.npred)]',[(oo_.dr.npred),1]);
%     %the foregoing takes into account that D is a kronecker and not block kronecker product, so to extract the kron([MA_.Y_i{j-1};zeros(M_.exo_nbr,M_.exo_nbr)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)]);
%     %for each j and i, well have to gather the scattered blocks.
    
    MA_.Y_ji{1,1}=(1/2)*MA_.BETA_00;
    for j=2:MA_.max_length
    MA_.Y_ji{j,1}=(1/2)*MA_.BETA_20*kron(MA_.Y_i{j-1}(select_state,:),eye(M_.exo_nbr));
    end 
    for i=2:MA_.max_length
    MA_.Y_ji{1,i}=MA_.Y_ji{i,1}*commute_ji_ij;%replacement code for: (1/2)*MA_.BETA_2*kron([zeros(M_.endo_nbr,M_.exo_nbr);eye(M_.exo_nbr)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)]);
     end
    for j=2:MA_.max_length
        for i=2:MA_.max_length
            %replacement code for: MA_.Y_ji{j,i}=MA_.ALPHA*MA_.Y_ji{j-1,i-1}+MA_.BETA_2*kron([MA_.Y_i{j-1};zeros(M_.exo_nbr,M_.exo_nbr)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)]);
            if i>=j
            MA_.Y_ji{j,i}=MA_.ALPHA*MA_.Y_ji{j-1,i-1}(select_state,:)+(1/2)*MA_.BETA_22*alt_kron(MA_.Y_i{j-1}(select_state,:),MA_.Y_i{i-1}(select_state,:));
            else
            MA_.Y_ji{j,i}=MA_.Y_ji{i,j}*commute_ji_ij;
            end    
        end
    end
    %Compute the third-order kernel
if options_.order == 3
    %preallocate
     K_matrix=commutation(M_.exo_nbr,M_.exo_nbr^2)*kron(eye(M_.exo_nbr),commutation(M_.exo_nbr,M_.exo_nbr));
     commute_kji_kij=kron(eye(M_.exo_nbr),commute_ji_ij);
     commute_kji_jki=kron(commute_ji_ij,eye(M_.exo_nbr));
     commute_kji_ikj=commutation(M_.exo_nbr,M_.exo_nbr^2);
     commute_kji_jik=commutation(M_.exo_nbr^2,M_.exo_nbr);
     [ia_ab,ib_ab]=meshgrid(1:M_.endo_nbr+M_.exo_nbr,1:(M_.endo_nbr+M_.exo_nbr)^2);
     [ja_ab,jb_ab]=meshgrid(1:M_.exo_nbr,1:(M_.exo_nbr)^2);
     [ic_ca,ia_ca]=meshgrid(1:M_.endo_nbr+(M_.endo_nbr+M_.exo_nbr)^2,1:(M_.endo_nbr+M_.exo_nbr));
     [jc_ca,ja_ca]=meshgrid(1:M_.exo_nbr^2,1:(M_.exo_nbr));
     [ia_ac,ic_ac]=meshgrid(1:(M_.endo_nbr+M_.exo_nbr),1:M_.endo_nbr+(M_.endo_nbr+M_.exo_nbr)^2);
     [ja_ac,jc_ac]=meshgrid(1:(M_.exo_nbr),1:M_.exo_nbr^2);
     
     MA_.Y_kji{1,1,1}=(1/6)*MA_.BETA_000;
     for k=2:MA_.max_length
        MA_.Y_kji{k,1,1}=(1/6)*MA_.BETA_300*alt_kron(MA_.Y_i{k-1}(select_state,:),eye(M_.exo_nbr^2));
       for i=2:MA_.max_length
       %an iterative replacement code for:   
       %MA_.Y_kji{k,1,i}=(1/6)*MA_.BETA_3*...
       %[kron([MA_.Y_i{k-1};zeros(M_.exo_nbr,M_.exo_nbr)],kron([zeros(M_.endo_nbr,M_.exo_nbr);eye(M_.exo_nbr)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)]));...
         %kron([MA_.Y_ji{k-1,1};kron([MA_.Y_i{k-1};zeros(M_.exo_nbr,M_.exo_nbr)],[zeros(M_.endo_nbr,M_.exo_nbr);eye(M_.exo_nbr)])],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)]);...
         %kron([zeros(M_.endo_nbr,M_.exo_nbr);eye(M_.exo_nbr)],[MA_.Y_ji{k-1,i-1};kron([MA_.Y_i{k-1};zeros(M_.exo_nbr,M_.exo_nbr)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)])])*K_matrix;...
         %kron([MA_.Y_i{k-1};zeros(M_.exo_nbr,M_.exo_nbr)],[MA_.Y_ji{1,i-1}
         %;kron([zeros(M_.endo_nbr,M_.exo_nbr);eye(M_.exo_nbr)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)])])];
        MA_.Y_kji{k,i,1}=(1/6)*MA_.BETA_330_1*alt_kron(alt_kron(MA_.Y_i{k-1}(select_state,:),MA_.Y_i{i-1}(select_state,:)),eye(M_.exo_nbr))+(1/6)*MA_.BETA_20*alt_kron(2*MA_.Y_ji{k,i}(select_state,:),eye(M_.exo_nbr));    
       end
         for j=2:MA_.max_length
             MA_.Y_kji{k,1,j}=MA_.Y_kji{k,j,1}*commute_kji_kij;
         end
     end
     for j=2:MA_.max_length
        MA_.Y_kji{1,j,1}=MA_.Y_kji{j,1,1}*commute_kji_jki;
      for i=2:MA_.max_length
        MA_.Y_kji{1,j,i}=MA_.Y_kji{j,1,i}*commute_kji_jki;
         end
     end
     for i=2:MA_.max_length
        MA_.Y_kji{1,1,i}=MA_.Y_kji{i,1,1}*commute_kji_ikj;
     end
     for k=2:MA_.max_length
         for j=2:MA_.max_length
             for i=2:MA_.max_length
                 if j<k
                    MA_.Y_kji{k,j,i}=MA_.Y_kji{j,k,i}*commute_kji_jki; 
                 elseif i<j
                     MA_.Y_kji{k,j,i}=MA_.Y_kji{k,i,j}*commute_kji_kij; 
                 else
                %an iterative replacement code for:  
                 % MA_.Y_kji{k,j,i}=MA_.ALPHA*MA_.Y_kji{k-1,j-1,i-1}+(1/6)*MA_.BETA_3*...
                 %   [repeated_kron([MA_.Y_i{k-1};zeros(M_.exo_nbr,M_.exo_nbr)],DkronD((M_.endo_nbr+M_.exo_nbr)*(MA_.max_length*(M_.endo_nbr+M_.exo_nbr))*(j-2)+(i-2)*(M_.endo_nbr+M_.exo_nbr)+select,:),ia_ab,ib_ab,ja_ab,jb_ab);...
                 %    repeated_kron([MA_.Y_ji{k-1,j-1};DkronD((M_.endo_nbr+M_.exo_nbr)*(MA_.max_length*(M_.endo_nbr+M_.exo_nbr))*(k-2)+(j-2)*(M_.endo_nbr+M_.exo_nbr)+select,:)],[MA_.Y_i{i-1};zeros(M_.exo_nbr,M_.exo_nbr)],ic_ca,ia_ca,jc_ca,ja_ca);...
                 %   repeated_kron([MA_.Y_i{j-1};zeros(M_.exo_nbr,M_.exo_nbr)],[MA_.Y_ji{k-1,i-1};DkronD((M_.endo_nbr+M_.exo_nbr)*(MA_.max_length*(M_.endo_nbr+M_.exo_nbr))*(k-2)+(i-2)*(M_.endo_nbr+M_.exo_nbr)+select,:)],ia_ac,ic_ac,ja_ac,jc_ac)*K_matrix;...
                  %  repeated_kron([MA_.Y_i{k-1};zeros(M_.exo_nbr,M_.exo_nbr)],[MA_.Y_ji{j-1,i-1};DkronD((M_.endo_nbr+M_.exo_nbr)*(MA_.max_length*(M_.endo_nbr+M_.exo_nbr))*(j-2)+(i-2)*(M_.endo_nbr+M_.exo_nbr)+select,:)],ia_ac,ic_ac,ja_ac,jc_ac)];
                        temp=alt_kron(2*MA_.Y_ji{k-1,j-1}(select_state,:),MA_.Y_i{i-1}(select_state,:))+alt_kron(2*MA_.Y_ji{k-1,i-1}(select_state,:),MA_.Y_i{j-1}(select_state,:))*commute_kji_kij+alt_kron(2*MA_.Y_ji{j-1,i-1}(select_state,:),MA_.Y_i{k-1}(select_state,:))*commute_kji_jik;
                      MA_.Y_kji{k,j,i}=MA_.ALPHA*MA_.Y_kji{k-1,j-1,i-1}(select_state,:)+(1/6)*MA_.BETA_333_1*alt_kron(MA_.Y_i{k-1}(select_state,:),alt_kron(MA_.Y_i{j-1}(select_state,:),MA_.Y_i{i-1}(select_state,:)))+(1/6)*MA_.BETA_22*temp;
                 end
             end
         end
     end
% Finally the correction to the first-order kernel
  MA_.Y_sigma_2_i{1}=(1/2)*MA_.BETA_sigma_2_0;                             
     for i=2:MA_.max_length
       MA_.Y_sigma_2_i{i}=MA_.ALPHA*MA_.Y_sigma_2_i{i-1}(select_state,:)+(1/2)*MA_.BETA_sigma_2_1*MA_.Y_i{i-1}(select_state,:);
         end
 end
end
