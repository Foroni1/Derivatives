clear all
clc

%% Import and manipulaiton of data (Option data)
	
% Import
data = readtable('20190102_SPY.csv');

% SPY price at 02/01/2019
S0 = data.UnderlyingPrice(1);
% T-Bill 3m at 02/01/2019 
r = 0.02417; 

%% Pre-processing

n = numel(categories(categorical(data.StrikePrice)));
m = numel(categories(categorical(data.ExpirationY)));

% Indexes of calls and puts
for i = 1:numel(data.PutCall)
    call_idx(i,1) = cellfun(@isequal,data.PutCall(i),{'call'});
    put_idx(i,1) = cellfun(@isequal,data.PutCall(i),{'put'});
end

% Extraction of strikes, maturities and prices
strikes = data.StrikePrice(call_idx);
mat = data.ExpirationY(call_idx);
call_prices = data.MidPrice(call_idx);
put_prices = data.MidPrice(put_idx);

strikes = reshape(strikes,[n m]);
mat = reshape(mat,[n m]);
call_prices = reshape(call_prices,[n m]);
put_prices = reshape(put_prices,[n m]);

money = strikes/S0*100;

%Plot price surface

% surf(strikes,mat,call_prices)
% surf(strikes,mat,put_prices)

%% Implied volatilities

call_imp_vol = blsimpv(S0*ones(n,m),strikes,r*ones(n,m),mat,call_prices);
call_imp_vol = fillmissing(call_imp_vol,'linear',2);

put_imp_vol = blsimpv(S0*ones(n,m),strikes,r*ones(n,m),mat,put_prices,'Class',{'put'});
put_imp_vol = fillmissing(put_imp_vol,'linear',2);

%% Volatility surface plot

%Call implied volatility surface
figure(1)
surf(money,mat,call_imp_vol)
title('SPY call volatility surface')
xlabel('Moneyness');
ylabel('Maturity (years)')
zlabel('Implied volatility')
colormap turbo
shading interp

%Put implied volatility surface
figure(2)
surf(money,mat,put_imp_vol,'EdgeColor','none')
title('SPY put volatility surface')
xlabel('Moneyness');
ylabel('Maturity (years)')
zlabel('Implied volatility')
colormap turbo
shading interp

%% Heston Calibration

% x0 = [kappa, theta, eta, rho, v0]
x0_c = [1 0.1 0.15 -0.01 0.03];
x0_p = [0.08 0.01 0.01 -0.01 0.02];

% Lower and upper bound for optimisation
lb = [0 0 0 -1 0];
ub = [20 1 5 0 1];

% Optimisation options
options = optimset('TolFun', 1e-8);

% Calibrations using non-linear least squares and a cost function
x_c = lsqnonlin(@(x) costf3(x,mat,strikes,S0*ones(n,m),call_imp_vol,money),x0_c,lb,ub,options);
x_p = lsqnonlin(@(x) costf3(x,mat,strikes,S0*ones(n,m),put_imp_vol,money),x0_p,lb,ub,options);

% Extraction of parameters
kappa_c = x_c(1); 
theta_c = x_c(2);
eta_c = x_c(3);
rho_c = x_c(4);
v0_c = x_c(5);

kappa_p = x_p(1); 
theta_p = x_p(2);
eta_p = x_p(3);
rho_p = x_p(4);
v0_p = x_p(5);

%% Heston simulation

% Options with maturity at 29/03/2019 (~3 months)
% Number of simulations
M = 10000;

% Time discretization
T = years(datetime(2019,3,29)-datetime(2019,1,2));
N = days(datetime(2019,3,29)-datetime(2019,1,2));
dt = T/N;

% Choleski factorization for calls and puts
A_c = chol([1 rho_c; rho_c 1]);
A_p = chol([1 rho_p; rho_p 1]);

% Random normal numbers
X_c = randn(M,N,2);
X_p = randn(M,N,2);

% Stock and volatility processes
Z_c = zeros(M,N+1);
V_c = zeros(M,N+1);
V_c(:,1) = v0_c;
Z_p = zeros(M,N+1);
V_p = zeros(M,N+1);
V_p(:,1) = v0_p;

for i = 2:N+1
    b_c = squeeze(X_c(:,i-1,:))*A_c;
    V_c(:,i) = abs(V_c(:,i-1)+kappa_c*(theta_c-V_c(:,i-1))*dt+eta_c*sqrt(V_c(:,i-1)*dt).*b_c(:,1));
    Z_c(:,i) = Z_c(:,i-1)+(r-V_c(:,i-1)/2)*dt+sqrt(V_c(:,i-1)*dt).*b_c(:,2);
    b_p = squeeze(X_p(:,i-1,:))*A_p;
    V_p(:,i) = abs(V_p(:,i-1)+kappa_p*(theta_p-V_p(:,i-1))*dt+eta_p*sqrt(V_p(:,i-1)*dt).*b_p(:,1));
    Z_p(:,i) = Z_p(:,i-1)+(r-V_p(:,i-1)/2)*dt+sqrt(V_p(:,i-1)*dt).*b_p(:,2);
end

S_c = S0*exp(Z_c);
S_p = S0*exp(Z_p);

%% Simulation Plot

t = linspace(0,T,N+1);

figure(3)
subplot(2,1,1);
plot(t,S_c);
title('Underlying simulation call');
xlabel('Time');
ylabel('Price');
subplot(2,1,2);
plot(t,V_c);
title('Volatility simulation call');
xlabel('Time');
ylabel('Volatility');

figure(4)
subplot(2,1,1);
plot(t,S_p);
title('Underlying simulation put');
xlabel('Time');
ylabel('Price');
subplot(2,1,2);
plot(t,V_p);
title('Volatility simulation put');
xlabel('Time');
ylabel('Volatility');

%% Pricing

% Import SPY time series
spy = readtable('SPY_data.xlsx');
adj_close = table2array(spy(:,2));

% Creation of returns and BS volatility
spy_ret = zeros(numel(adj_close),1);
for i = 2:numel(adj_close)
    spy_ret(i,1) = log(adj_close(i)/adj_close(i-1));
end

sig_bs = std(spy_ret)*sqrt(250);

% Payoffs calls and puts with Heston
K = 250;
payoff_c = max(S_c(:,end)-K,0);
payoff_p = max(K-S_p(:,end),0);

% Options prices with Heston and BS

put_price_heston = mean(payoff_p);
call_price_heston = mean(payoff_c);

call_price_bs = callBS(S0,K,r,sig_bs,T);
put_price_bs = putBS(S0,K,r,sig_bs,T);

% Print prices
fprintf("Heston's model call price is: %3.2f\n",call_price_heston)
fprintf("Heston's model put price is: %3.2f\n",put_price_heston)
fprintf("BS's model call price is: %3.2f\n",call_price_bs)
fprintf("BS's model put price is: %3.2f\n",put_price_bs)
