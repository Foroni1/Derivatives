function P=putBS(S,K,r,sig,T)
% Calculates the price of an European put option under the Black and Scholes model

% INPUT: 
%   S   = price of the underlying
%   K   = strike
%   r   = interest rate
%   sig = volatility
%   T   = expiration

% OUTPUT:
%   P   = price of the put option

d1 = (log(S/K)+(r+sig^2/2)*T)/(sig*sqrt(T));
d2 = d1-sig*sqrt(T);
P  = K*exp(-r*T)*normcdf(-d2)-S*normcdf(-d1);
