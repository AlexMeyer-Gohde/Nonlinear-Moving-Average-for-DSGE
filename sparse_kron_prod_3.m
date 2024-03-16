function [FF]=sparse_kron_prod_3(DD,CC,AA,BB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sparse_kron_prod_3.m
%
% This file calculates DD*kron(C,A,B) by only performing nonzero operation. 
% Currently, there is no mex equivalent from Dynare for this operation.
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
[uu ooo]=size(CC);
[tt nnn]=size(AA);
[ss mmm]=size(BB);
   

[df_row,df_col]=size(DD);
[row,col,DD_nonzero] = find(DD);
[CC_row,CC_col,CC_nonzero] = find(CC);
[AA_row,AA_col,AA_nonzero] = find(AA);
[BB_row,BB_col,BB_nonzero] = find(BB);
[row,rowi] = sort(row);col=col(rowi);DD_nonzero=DD_nonzero(rowi);

FF=zeros(df_row,ooo*nnn*mmm);
A_kron_B=zeros(1,nnn*mmm);
C_kron_A_kron_B=zeros(1,ooo*nnn*mmm);
FFF=zeros(1,ooo*nnn*mmm);
for i=1:df_row
    FFF=zeros(1,ooo*nnn*mmm);
    col_op=col(row==i);
    if isempty(col_op)==0
         col_mult=DD_nonzero(row==i);
                    ll=ceil(col_op./(ss*tt));
           ll_2=ceil(col_op./ss-(ll-1).*tt);
           kk=col_op-(ll-1).*ss*tt-(ll_2-1).*ss;
           
        for j=1:length(col_op)     
    oo=CC_col(CC_row==ll(j));
    nn=AA_col(AA_row==ll_2(j));
    mm=BB_col(BB_row==kk(j));
            if isempty(oo)==0&&isempty(nn)==0&&isempty(mm)==0
      l_oo=length(oo);l_nn=length(nn);l_mm=length(mm); 
      oo1=oo(:,ones(1,l_nn*l_mm)).';oo1=oo1(:);
      nn1=nn(:,ones(1,l_mm)).';nn1=nn1(:);inn1=(1:l_mm*l_nn).'; inn1 = inn1(:,ones(1,l_oo));nn1=nn1(inn1,1);%nn1=repmat(nn1,[l_oo 1]);
      imm1=(1:l_mm).'; imm1 = imm1(:,ones(1,l_nn*l_oo));mm1=mm(imm1,1);%mm1=repmat(mm,[l_nn*l_oo 1]);
      bb=mm1+(oo1-1).*(nnn*mmm)+(nn1-1).*mmm;   
A_kron_B=zeros(1,l_nn*l_mm);
C_kron_A_kron_B=zeros(1,l_oo*l_nn*l_mm);
      A_kron_B=alt_kron(AA(ll_2(j),nn),BB(kk(j),mm));
      C_kron_A_kron_B=alt_kron(CC(ll(j),oo),A_kron_B);
           FFF(bb)=FFF(bb)+col_mult(j)*C_kron_A_kron_B;
            end
        end
        FF(i,:)=FFF;
    end
end