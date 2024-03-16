function [MA_]=first_order(M_,MA_,oo_,options_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% first_order.m
%
% This file produces the coefficients for a first-order approximation 
% 
% solution form: b(i) = ALPHA * b(i-1) + BETA_1ST * U(i), where U(i) = 0^i
%
% If only a first or second order approximation was calculated, use dynares
% decision rules
%
% If a third order approximation was calculated by dynare, we must resolve
% the first order to extract the certainty equivalent component of the
% first derivative of the policy function.
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
% Copyright (C) 2005-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if options_.order>=3
    % Recover certainty equivalent dr
    options_.order=1;
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;%not always set by dynare for higher order perturbations
    end
    exo_simul=oo_.exo_simul;oo_.exo_simul=zeros(M_.maximum_lag+2,M_.exo_nbr);% these two lines are needed to prevent dynare from using the simulation to evaluate the derivatives
    exo_det_simul=oo_.exo_det_simul;exo_det_simul=[];
    if MA_.dynare_version<4.3 %Starting in Dynare 4.3.0, dr1 is removed. 
        [oo_.dr,waste1,waste2,waste3,waste4] = dr1(oo_.dr,0,M_,options_,oo_);
    elseif MA_.dynare_version>=4.3 %replace with stochastic_solvers in 4.3
        [oo_.dr,waste1] =stochastic_solvers(oo_.dr,0,M_,options_,oo_);
    else
        disp('Error, no certainty evuivalent dr calculated');
    end
    options_.order=3;oo_.exo_simul=exo_simul;oo_.exo_det_simul=exo_det_simul;
end
MA_.ALPHA=oo_.dr.ghx;
MA_.BETA_0=oo_.dr.ghu;




end

