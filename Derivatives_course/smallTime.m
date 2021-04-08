function H = SmallTime(kappa,theta,sigma,rho,v0,money,T)
%This function is used within the Heston calibration routine to speed up the process 
%using an asymptotic formula. Precisely, this function has to be used in a cost function
%during the calibration

%INPUT:
%   kappa: scalar rate at which the variance returns to its long run mean
%   theta: scalar long run mean of the variance
%   sigma: scalar volatility of the variance
%   rho: scalar correlation between the price process and the variance process
%   v0: scalar initial variance
%   money: scalar/matrix of options moneyness
%   T: scalar/matrix of options maturities

%OUTPUT:
%   H: scalar/matrix of options asymptotical implied volatilities under the Heston model

x = log(money);
alpha = kappa*theta;
term1 = 1 + (rho/2)*(sigma*x/v0) + (1-(7/4)*rho^2)* (sigma^2*x^2)/(12*v0^2);
term2a = (rho*sigma*v0)/2 - (sigma^2/6)*(1-rho^2/4) + alpha - kappa*v0;
term2b = sigma^2*(1-rho^2) + rho*sigma*v0 - 2*alpha - 2*kappa*v0;
term3a = (176 - 712*rho^2 + 521*rho^4)*sigma^2 + 40*sigma*rho^3*v0;
term3b = 80*(13*rho^2 - 6)*alpha - 80*kappa*rho^2*v0;
H1 = v0 * term1 + (term2a + (rho*sigma)/(12*v0) * term2b * x) * (T/2);
H2 = (sigma^2/(7680*v0^2)) * (term3a + term3b) * x^2*T;
H = H1 + H2;
end
