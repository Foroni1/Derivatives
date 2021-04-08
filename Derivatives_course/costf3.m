function [L] = costf3(x, T, strike, Under, impl_vol, money)

%Cost function used in the Heston calibration

%INPUT:
%   x: vector of initial parameters
%   T: scalar/vector of maturities
%   strike: scalar/vector of options strikes
%   Under: scalar/matrix of underlying prices (must have the same rows of strike and the same columns of T)
%   impl_vol: scalar/matrix of options implied volatilities (same size of Under)
%   money: scalar/vector of options moneyness

%OUTPUT:
%   L: residual value of the cost function

for i = 1 : size(strike,1)
    for j = 1:size(strike,2)
        L(i,j)=((impl_vol(i,j)^2-SmallTime(x(1),x(2),x(3),x(4),x(5),money(i,j),T(i,j)))/impl_vol(i,j)^2);
    end
end
end
