%Model of Aruoba et al 2006 JEDC

var k c l z;
varexo e;

parameters alpha beta delta theta tau rho sigma;

alpha  = 0.4;
beta   = 0.9896; 
delta  = 0.0196; 
theta  = 0.357; 
tau    = 10; 
rho    = 0.95; 
sigma  = 0.007;

model; 
  (((exp(c)^theta)*(1-exp(l))^(1-theta))^(1-tau))/exp(c) = beta*((((exp(c(+1))^theta)*(1-exp(l(+1)))^(1-theta))^(1-tau))/exp(c(+1)))*(1-delta+alpha*exp(z(+1))*(exp(k)^(alpha-1))*exp(l(+1))^(1-alpha));
exp(c)+exp(k)= exp(z)*(exp(k(-1))^(alpha))*exp(l)^(1-alpha)+(1-delta)*exp(k(-1));
  (((exp(c)^theta)*(1-exp(l))^(1-theta))^(1-tau))*(1-theta)/(1-exp(l))=(((exp(c)^theta)*(1-exp(l))^(1-theta))^(1-tau))*theta*(1-alpha)*exp(z)*exp(k(-1))^alpha*exp(l)^(-alpha)/exp(c);
  z = rho*z(-1)+e;
end;

steady_state_model;
phi     = (1/alpha*(1/beta-1+delta))^(1/(1-alpha));
omega   = (phi^(1-alpha)-delta);
psi = theta/(1-theta)*(1-alpha)*phi^(-alpha);
  
k = log(psi/(omega+phi*psi));
l = log(phi*exp(k));
c = log(omega*exp(k));
z=0;
end;

steady;

shocks;
var e = sigma^2;
end;

stoch_simul(irf=0, order = 3, periods=200, drop=50);
options_.irf=100;
my_nonlinear_MA_options.plot_simulations=0;
my_nonlinear_MA_options.shock_scale=-10;
my_nonlinear_MA_options.plot_irf=0;
[MA_]=nonlinear_MA(M_,oo_,options_,my_nonlinear_MA_options);
