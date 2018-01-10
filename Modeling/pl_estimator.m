function [ pl_db ] = pl_estimator(param,d_min,d_max )
%PL_ESTIMATOR Summary of this function goes here
%   Detailed explanation goes here
distance= d_min:d_max;
Tx_height = param.tx_h;
Rx_height = Tx_height;
alpha = param.alpha;
Epsilon_r = param.eps;
CarrierFreq=param.freq;
Speedof_Light=3*10^8;



lambda=Speedof_Light/CarrierFreq;
Path_Loss_dB = zeros(d_max,1);
Mean_of_Fading = zeros(d_max,1);
dbmean = zeros(d_max,1);
m_meanr2 = zeros(d_max,1);
m_meanr3 = zeros(d_max,1);
k_hat = zeros(d_max,1);
mu_hat = zeros(d_max,1);
alpha_hat = zeros(d_max,1);
omega_hat = zeros(d_max,1);
m_hat = zeros(d_max,1);
etha_hat = zeros(d_max,1);
u3_hat = zeros(d_max,1);
sigma2_h_Fat = zeros(d_max,1);
NakMu_hat = zeros(d_max,1);
NakOmega_hat = zeros(d_max,1);

distance_los = sqrt(distance.^2+(Tx_height-Rx_height)^2);
distance_ref = sqrt(distance.^2+(Tx_height+Rx_height)^2);

sint = (Tx_height+Rx_height)./distance_ref;
cost = distance./distance_ref;

Gamma = (sint-sqrt(Epsilon_r-cost.^2))./(sint+sqrt(Epsilon_r-cost.^2));
Phi = 2*pi*(distance_los-distance_ref)/lambda;

temp = 1+(Gamma .* exp(1i*Phi));
Path_Loss_dB(d_min:d_max) = 10*alpha*log10((4*pi)*(distance/lambda).*(abs(temp).^(-1)));
pl_db = Path_Loss_dB;
end

