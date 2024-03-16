 function [MA_]=plots_kernels(M_,MA_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plots_kernels.m
%
% This file produces plots of the Volterra kernels calculated in 
% kernels.m
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
 if options_.order >= 2
    temp = cat(3,MA_.Y_ji{:,:});
    for jj=1:M_.exo_nbr
        for kk=1:jj
            temp_2 = reshape(squeeze(temp(:,(jj-1)*M_.exo_nbr+kk,:,:)), [M_.endo_nbr, MA_.max_length, MA_.max_length]);
            for ii=1:M_.endo_nbr
                if max(max(squeeze(temp_2(ii,:,:))))~=min(min(squeeze(temp_2(ii,:,:))))
                figure;
                 mesh(0:MA_.max_length-1,0:MA_.max_length-1,squeeze(temp_2(ii,:,:)));
                eval(sprintf('title(''Second-Order Kernel of %s with Respect to %s and %s'')',M_.endo_names(oo_.dr.order_var(ii),:),M_.exo_names(jj,:),M_.exo_names(kk,:)))
                eval(sprintf('zlabel(''Second-Order Component of %s at Time t'')',M_.endo_names(oo_.dr.order_var(ii),:)))
                eval(sprintf('xlabel(''Shock in %s at Time t-i'')',M_.exo_names(kk,:)))
                eval(sprintf('ylabel(''Shock in %s at Time t-j'')',M_.exo_names(jj,:)))
                end
            end
        end
    end
     if options_.order == 3
         temp = cat(3,MA_.Y_kji{:,:,:});
            for jj=1:M_.exo_nbr
                temp_2 = reshape(squeeze(temp(:,(jj-1)*M_.exo_nbr^2+(jj-1)*M_.exo_nbr+jj,:,:)), [M_.endo_nbr, MA_.max_length, MA_.max_length, MA_.max_length]);
                for ii=1:M_.endo_nbr
                    if MA_.plot_kernels==1
                       temp_3=squeeze(temp_2(ii,:,:,:));
                        if max(max(max(temp_3)))~=min(min(min(temp_3)))
                             figure;colormap(jet(36))
                        [x,y,z] = meshgrid(0:1:MA_.max_length-1,0:1:MA_.max_length-1,0:1:MA_.max_length-1);
                        hsp = surf(linspace(0,MA_.max_length-1,MA_.max_length),linspace(0,MA_.max_length-1,MA_.max_length),0.5*(MA_.max_length)*ones(MA_.max_length));
                        rotate(hsp,[1,-1,0],45);
                         xd = get(hsp,'XData');
                         yd = get(hsp,'YData');
                         zd = get(hsp,'ZData');
                         delete(hsp)
                            h=slice(x,y,z,temp_3,xd,yd,zd);
                            colorbar
                            set(h,'EdgeColor','none')
                            eval(sprintf('title(''Third-Order Kernel of %s with Respect to %s'')',M_.endo_names(oo_.dr.order_var(ii),:),M_.exo_names(jj,:)))
                            eval(sprintf('xlabel(''Shock in %s at Time t-i'')',M_.exo_names(jj,:)))
                            eval(sprintf('ylabel(''Shock in %s at Time t-j'')',M_.exo_names(jj,:)))
                            eval(sprintf('zlabel(''Shock in %s at Time t-k'')',M_.exo_names(jj,:)))
                            hold on
                            caxis auto
                            color_lim = caxis;
                            cont_intervals = linspace(color_lim(1),color_lim(2),10);
                            hcont=contourslice(x,y,z,temp_3,xd,yd,zd,cont_intervals,'linear');
                            set(hcont,'EdgeColor',[.4 .4 .4],'LineWidth',1)
                            hx=slice(x,y,z,temp_3,MA_.max_length-1,[],[]);
                            set(hx,'EdgeColor','none')
                            hy=slice(x,y,z,temp_3,[],MA_.max_length-1,[]);
                            set(hy,'EdgeColor','none')
                            hz=slice(x,y,z,temp_3,[],[],0);
                            set(hz,'EdgeColor','none')
                            contourslice(x,y,z,temp_3,MA_.max_length-1,[],[]);
                            contourslice(x,y,z,temp_3,[],MA_.max_length-1,[]);
                            contourslice(x,y,z,temp_3,[],[],0);
                            hold off
                         end
                    end
                end
            end
        end
 end