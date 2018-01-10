clc
clearvars -except data_set
close all
%% preparation

class param
input  = 'same_density_0to12.csv';
if ~exist('data_set','var')==1
    csv_data = readtable(input);
    data_set = [csv_data.Range,csv_data.RSS];
    data_set=double(data_set);
    [r,IX]=sort(data_set(:,1));
    data_set(:,1)=r;
    data_set(:,2)=data_set(IX,2);
end
% plot(data_set(:,1),data_set(:,2),'.')
%% Dataset Mean of RSSI-Distance Plot

d_max = 89;
d_min = 10;

fixed_fading_distance = 10;

Tx_Power = 20;

index=[];

RSSI_per_meter_mean = zeros(d_max,1);
total_point_per_meter = zeros(d_max,1);
% total_p = zeros(d_max-d_min+1,1);
lower = find(data_set(:,1)<=d_min-1,1,'last');

for k=d_min:d_max
    
    upper = find(data_set(:,1)<=k,1,'last');
    RSSI_per_meter_mean(k,1)=mean(data_set(lower+1:upper,2));
    if( ~isnan(RSSI_per_meter_mean(k,1)))
        index = [index,k];
    end
    %     ms(i,2)=std(data_set(lower+1:upper,2));
    total_point_per_meter(k) = upper - lower;
    lower = upper;
end

clear lower upper

% plot(index,RSSI_per_meter_mean(index),'.')

%% Parameters
Speedof_Light=3*10^8;
alpha=2.02
CarrierFreq=5.2649*10^9;
Tx_height = 1.565
Epsilon_r = 1.0148
Rx_height = Tx_height;
lambda=Speedof_Light/CarrierFreq;
%% Calc the path-loss

distance= d_min:d_max;

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
% Path loss calculated
lower = find(data_set(:,1)<=d_min-1,1,'last');
%% Calculating Linear Fading of Dataset linear_fds=(sqrt(10^fading)/10)/1000 and Linear = (sqrt(10^RSSI)/10)/1000
for k=d_min:d_max
    
    %% linear and dBm fading only
    upper = find(data_set(:,1)<=k,1,'last');
    fds(lower+1:upper,2) = data_set(lower+1:upper,2)-Tx_Power+Path_Loss_dB(k);
    
    fds(lower+1:upper,1) = data_set(lower+1:upper,1);
    
   linear_fds(lower+1:upper,2) = sqrt((10.^((fds(lower+1:upper,2)) ./ 10))./1000);
    
   
   linear_fds(lower+1:upper,1) = fds(lower+1:upper,1);
    
   linear_data(lower+1:upper,2) = sqrt((10.^((data_set(lower+1:upper,2)) ./ 10))./1000);
   
   linear_data(lower+1:upper,1) = fds(lower+1:upper,1);
    lower = upper;
end
clear lower upper


lower = find(data_set(:,1)<=d_min-1,1,'last');

for k=d_min+fixed_fading_distance-1:fixed_fading_distance:d_max
    
    upper = find(data_set(:,1)<=k,1,'last');
    
    %% per meter m calculations
    dbmean(k) = mean(fds(lower+1:upper,2));
    m_red_fds(lower+1:upper,1) = fds(lower+1:upper,2) - dbmean(k);
    m_red_2(lower+1:upper,1) = m_red_fds(lower+1:upper).^2;
    m_red_3(lower+1:upper,1) = m_red_fds(lower+1:upper).^3;
    m_meanr2(k) = mean(m_red_2(lower+1:upper,1));
    m_meanr3(k) = mean(m_red_3(lower+1:upper,1));
    k_hat(k) = m_meanr2(k)^(3.0/2.0)/m_meanr3(k);
    
    for j=1:fixed_fading_distance-1
        k_hat(k-j) = k_hat(k);
    end
    
    %% Alpha-Mu
    if ( k_hat(k) <= -2.85)
        mu_hat(k) = k_hat(k)^2+0.5;
        alpha_hat(k) = (20/log(10))*sqrt(psi(1,mu_hat(k))/m_meanr2(k));
        omega_hat(k) = (mean(linear_fds(lower+1:upper,2).^alpha_hat(k))).^(1/alpha_hat(k));
        
        for j=1:fixed_fading_distance-1
            k_hat(k-j) = k_hat(k);
            mu_hat(k-j) = mu_hat(k);
            alpha_hat(k-j) = alpha_hat(k);
            omega_hat(k-j) = omega_hat(k);
        end
        
    elseif ((-2.85 < k_hat(k)) && (k_hat(k) <= -0.6))
        mu_hat(k) = -0.0773*k_hat(k)^4-0.6046*k_hat(k)^3-0.7949*k_hat(k)^2-2.4675*k_hat(k)-0.9208;
        alpha_hat(k) = (20/log(10))*sqrt(psi(1,mu_hat(k))/m_meanr2(k));
        omega_hat(k) = (mean(linear_fds(lower+1:upper,2).^alpha_hat(k))).^(1/alpha_hat(k));
        
        for j=1:fixed_fading_distance-1
            k_hat(k-j) = k_hat(k);
            mu_hat(k-j) = mu_hat(k);
            alpha_hat(k-j) = alpha_hat(k);
            omega_hat(k-j) = omega_hat(k);
        end
    elseif (-0.6 < (k_hat(k))&&( k_hat(k) < -0.5))
        mu_hat(k) = -132.8995*k_hat(k)^3-232.0659*k_hat(k)^2-137.6303*k_hat(k)-27.3616;
        alpha_hat(k) = (20/log(10))*sqrt(psi(1,mu_hat(k))/m_meanr2(k));
        omega_hat(k) = (mean(linear_fds(lower+1:upper,2).^alpha_hat(k))).^(1/alpha_hat(k));
        
        for j=1:fixed_fading_distance-1
            k_hat(k-j) = k_hat(k);
            mu_hat(k-j) = mu_hat(k);
            alpha_hat(k-j) = alpha_hat(k);
            omega_hat(k-j) = omega_hat(k);
        end
    else
        sprintf('invalid value for k_hat %f in bin %d',k_hat(k),k)
        mu_hat(k) = mu_hat(k-fixed_fading_distance);
        alpha_hat(k) = alpha_hat(k-fixed_fading_distance);
        omega_hat(k) = omega_hat(k-fixed_fading_distance);
        
        for j=1:fixed_fading_distance-1
            k_hat(k-j) = k_hat(k);
            mu_hat(k-j) = mu_hat(k);
            alpha_hat(k-j) = alpha_hat(k);
            omega_hat(k-j) = omega_hat(k);
        end
    end
    
    
    %Nakagami Distribution
    
    phat = mle(linear_fds(lower+1:upper,2),'distribution','nakagami');
    NakMu_hat(k) = phat(1);
    NakOmega_hat(k) = phat(2);
    
    for j=1:fixed_fading_distance-1
        NakMu_hat(k-j) = NakMu_hat(k);
        NakOmega_hat(k-j) = NakOmega_hat(k);
    end
    
    
    lower = upper;
    
end

clear lower upper

%% Fixing first bin
% k=13;
% mu_hat(k) = mu_hat(k+1);
% alpha_hat(k) = alpha_hat(k+1);
% omega_hat(k) = omega_hat(k+1);
% 
% for j=1:fixed_fading_distance-1
%     k_hat(k-j) = k_hat(k);
%     mu_hat(k-j) = mu_hat(k);
%     alpha_hat(k-j) = alpha_hat(k);
%     omega_hat(k-j) = omega_hat(k);
% end

%% alpha-mu and Nakagami random generator
rands = zeros(sum(total_point_per_meter),2);
repeat = 1;
lim = 0;
lowerf = find(linear_fds(:,1)<=d_min-1,1,'last');
hval = 5*ones(length(1:d_max),repeat);
pval = ones(length(1:d_max),repeat);
ksstat = ones(length(1:d_max),repeat);
Nakhval = 5*ones(length(1:d_max),repeat);
Nakpval = ones(length(1:d_max),repeat);
Nakksstat = ones(length(1:d_max),repeat);
alpha_mu_samples_res =[]
nakagami_samples_res = []
for z = d_min:d_max
    
    alpha_mu_samples = zeros(total_point_per_meter(z),repeat);
    Nakrands_ln = zeros(total_point_per_meter(z),repeat);
    sim = zeros(total_point_per_meter(z),repeat);
    Naksim = zeros(total_point_per_meter(z),repeat);
    
    upperf = find(linear_fds(:,1)<=z,1,'last');
    rands(lim+1:lim+total_point_per_meter(z),1) = z;
    
    %per meter
    a  = mu_hat(z);
    c = alpha_hat(z);
    b = omega_hat(z)/(a^(1/c));
    
    for s = 1:repeat
        alpha_mu_samples(:,s) = wggamrnd(a,b,c,total_point_per_meter(z),1);
    end
    e = NakMu_hat(z);
    f = NakOmega_hat(z);
    
    for s = 1:repeat
        Nakrands_ln(:,s) =  random('nakagami',e,f,total_point_per_meter(z),1);
    end
    
    %% KSTEST per meter for alpha-mu
        
    field = linear_fds(lowerf+1:upperf,2);
    for s=1:repeat
       sim(:,s) = alpha_mu_samples(:,s);
       Naksim(:,s) = Nakrands_ln(:,s);
       [hval(z,s),pval(z,s),ksstat(z,s)] = kstest2(field,sim(:,s),'Alpha',0.001);
       [Nakhval(z,s),Nakpval(z,s),Nakksstat(z,s)] = kstest2(field,Naksim(:,s),'Alpha',0.001);

    end
    d_temp = z*ones(total_point_per_meter(z),1)
    alpha_mu_samples_res = [alpha_mu_samples_res;[d_temp,alpha_mu_samples(:)]];
    nakagami_samples_res = [nakagami_samples_res;[d_temp,Nakrands_ln(:)];
  %  end
    %%
    lim = lim+total_point_per_meter(z);
    lowerf = upperf;
end



%% Adding path loss component

rnds = rands(rands(:,1)~=0,:);
lower = 0;

for k=d_min:d_max
    
    upper = find(rnds(:,1)<=k,1,'last');
    recons(lower+1:upper,2) = rnds(lower+1:upper,2)+Tx_Power-Path_Loss_dB(k);
    
    Mean_of_Fading= mean(rnds(lower+1:upper,2));
    
    recons(lower+1:upper,1) = rnds(lower+1:upper,1);
    
    lower = upper;
end

clear lower upper

tmpindx = find(recons(:,2) < -94);
if(~isempty(tmpindx))
    recons(tmpindx,2) = 0;
    recons(tmpindx,1) = 0;
end

%     recons(isnan(recons(:,2)),1) = 0;

recons = recons(recons(:,1)~=0,:);


Recv_Power = Tx_Power-Path_Loss_dB + Mean_of_Fading;

% figure
% plot(recons(:,1),recons(:,2),'b.',data_set(:,1),data_set(:,2),'r.')
% title('reconstruction')
% legend('reconstruction','field')
% hold on
% plot(data_set(:,1),data_set(:,2),'r.')


%% median plotting

figure

lower = find(data_set(:,1)<=d_min-1,1,'last');
upper = find(data_set(:,1)<=d_max,1,'last');

ax1 = boxplot(data_set(lower+1:upper,2),data_set(lower+1:upper,1),'plotstyle','compact','colors','r','outliersize',10^-10);

hold on
boxplot(round(recons(:,2)),recons(:,1),'plotstyle','compact','colors','b','outliersize',10^-10,'colors','b');

% fig = gcf;
% ax = fig.CurrentAxes;
% 
% set(ax,'XTickLabel',{' '})
% 
% h = get(ax,'children');
% 
% for j=1:length(h)
%     child = get(h(j),'Children');
%     texts = findobj(child,'Type','Text');
%     for i=1:length(texts)
%         texts(i).String = '';
%     end
% 
% end
% 
% xData = linspace(0,400,17);% hardcoded values
% ax.XTick = xData;
% ax.XTickLabel = xData;
% 
% hval(:,2)=1:d_max;
% ind1 = find(hval(:,1)==1);
% 
% 
% plot(hval(ind1(d_min:end),2),-20+(hval(ind1(d_min:end),1)-1),'m.')
% ind2 = find(hval(:,1)==0);
% 
% plot(hval(ind2(d_min:end),2),-90+(hval(ind2(d_min:end),1)),'k.')

%%

% 
% color = ['r', 'b'];
% 
% tb = findobj(ax,'Tag','Box');
% 
% 
% for j=1:2
%    patch(get(tb(j),'XData'),get(tb(j),'YData'),color(j),'FaceAlpha',.5);
% end
% 
% 
% 
% h = get(ax, 'Children');
% 
% hleg1 = legend(h(1:2), 'Feature1', 'Feature2' );













































































