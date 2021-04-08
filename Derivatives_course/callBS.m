function C=callBS(S,K,r,sig,T)
% Function that calculates the value of an European call option under the Black and Scholes model
  
% INPUT: 
%   S   = price of the underlying
%   K   = call strike price
%   r   = interest rate
%   sig = volatility 
%   T   = expiration
  
% OUTPUT:
%   C   = price of the call option

d1 = (log(S/K)+(r+sig^2/2)*T)/(sig*sqrt(T));
d2 = d1-sig*sqrt(T);
C  = S*normcdf(d1)-K*exp(-r*T)*normcdf(d2);
