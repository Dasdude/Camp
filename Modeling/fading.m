clc
clearvars -except data_set
close all
%% preparation
clear
input  = 'same_density_30to50.csv';
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

d_max = 700;
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
alpha=2.02;
CarrierFreq=5.89*10^9;
Tx_height = 1.565;
Epsilon_r = 1.0148;


%% params for freq 5.89 density 0-12
alpha=2.0108
CarrierFreq=5.89*10^9;
Tx_height = 1.4787
Epsilon_r = 1.0043
%% params for freq 5.89 density 12-130
alpha=2.0237
CarrierFreq=5.89*10^9;
Tx_height = 1.4787
Epsilon_r = 1.0091
%% params for freq 5.89 density 30-50
alpha=2.0576
CarrierFreq=5.89*10^9;
Tx_height = 1.4787
Epsilon_r = 1.0036
%% Fixed 
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
% Path_Loss_dB = RSSI_per_meter_mean
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
alpha_mu_samples_res =[];
nakagami_samples_res = [];
alpha_mu_samples_trunc =[];
nakagami_samples_trunc = [];
kld_alpha = zeros(size(1:d_max));
kld_nak = zeros(size(1:d_max));
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
        if total_point_per_meter(z)==0
            alpha_mu_samples(:,s)=ones(0,1)
            continue
        end
        alpha_mu_samples(:,s) = wggamrnd(a,b,c,total_point_per_meter(z),1);
    
    end
    e = NakMu_hat(z);
    f = NakOmega_hat(z);
    
    for s = 1:repeat
        Nakrands_ln(:,s) =  random('nakagami',e,f,total_point_per_meter(z),1);
    end
    
    %% KSTEST per meter for alpha-mu
    z
    field = linear_fds(lowerf+1:upperf,2);
    for s=1:repeat
       sim(:,s) = alpha_mu_samples(:,s);
       if all(isnan(sim(:,s)))
           continue
       end
       Naksim(:,s) = Nakrands_ln(:,s);
       [hval(z,s),pval(z,s),ksstat(z,s)] = kstest2(field,sim(:,s),'Alpha',0.01);
       [Nakhval(z,s),Nakpval(z,s),Nakksstat(z,s)] = kstest2(field,Naksim(:,s),'Alpha',0.001);
       %% KLD 
       binEdges    =  [-inf ; sort(field(:,s)) ; inf];
       nak_trunc_index = (10*log10((Naksim(:,s).^2)*1000)+Tx_Power-Path_Loss_dB(z)>=-94);
       am_trunc_index = (10*log10((sim(:,s).^2)*1000)+Tx_Power-Path_Loss_dB(z)>=-94);
       nak_trunc = Naksim(nak_trunc_index);
%        min(10*log10((nak_trunc(:,s).^2)*1000)+Tx_Power-Path_Loss_dB(z))
       am_trunc = sim(am_trunc_index);
       if size(am_trunc,1)==0
           am_trunc=ones(0,1)
       end
       if size(nak_trunc,1)==0
           nak_trunc=ones(0,1)
       end
       
       fieldc = histc(field,binEdges);
       nakc = histc(nak_trunc(:,s),binEdges);
       if size(nakc,2)>1
           nakc=nakc'
       end
       alphac = histc(am_trunc(:,s),binEdges);
       if size(alphac,2)>1
           alphac=alphac'
       end
       if size(fieldc,2)>1
           fieldc=fieldc'
       end
       nak_pdf = nakc./sum(nakc);
       alpha_pdf = alphac./sum(alphac);
       field_pdf = fieldc./sum(fieldc);
       bitdif_alpha = field_pdf.*log2((field_pdf./(alpha_pdf+eps))+eps);
       bitdif_nak = field_pdf.*log2((field_pdf./(nak_pdf+eps))+eps);
       kld_alpha(z) = sum(bitdif_alpha);
       kld_nak(z) = sum(bitdif_nak);
    end
    d_temp = z*ones(total_point_per_meter(z),1);
    alpha_mu_samples_res = [alpha_mu_samples_res;[d_temp,alpha_mu_samples(:)]];
    nakagami_samples_res = [nakagami_samples_res;[d_temp,Nakrands_ln(:)]];
%     alpha_mu_samples_trunc = [alpha_mu_samples_trunc;[z*ones(size(am_trunc(:),1),1),am_trunc(:)]];
%     nakagami_samples_trunc = [nakagami_samples_trunc;[z*ones(size(nak_trunc(:),1),1),nak_trunc]];
  %  end
    %%
    lim = lim+total_point_per_meter(z);
    lowerf = upperf;
end



%% Adding path loss component

% rnds = rands(rands(:,1)~=0,:);
lower = 0;
%% ALPHA MU RSSI Reconstruction
alpha_mu_recons = zeros(size(alpha_mu_samples));

alpha_mu_dbm = [alpha_mu_samples_res(:,1),10*log10((alpha_mu_samples_res(:,2).^2)*1000)];
alpha_mu_fading_mean= zeros(size(Path_Loss_dB));
for k=d_min:d_max
    
    upper = find(alpha_mu_samples_res(:,1)<=k,1,'last');
    alpha_mu_fading_mean(k)= mean(alpha_mu_dbm(lower+1:upper,2));
    alpha_mu_recons(lower+1:upper,2) = alpha_mu_dbm(lower+1:upper,2)+Tx_Power-Path_Loss_dB(k);
    alpha_mu_recons(lower+1:upper,1) = k;
    alpha_mu_fading_mean(k)= mean(alpha_mu_dbm(lower+1:upper,2));
    
    lower = upper;
end
%% NAKAGAMI RSSI Reconstruction
nakagami_recons = zeros(size(nakagami_samples_res));
nakagami_dbm = [nakagami_samples_res(:,1),10*log10((nakagami_samples_res(:,2).^2)*1000)];
nakagami_fading_mean = zeros(size(Path_Loss_dB));
for k=d_min:d_max
    upper = find(nakagami_samples_res(:,1)<=k,1,'last');
    nakagami_recons(lower+1:upper,2) = nakagami_dbm(lower+1:upper,2)+Tx_Power-Path_Loss_dB(k);
    nakagami_recons(lower+1:upper,1) = k;
    nakagami_fading_mean(k)= mean(nakagami_dbm(lower+1:upper,2));
    
    lower = upper;
end
clear lower upper


rx_power_nakagami_mean = Tx_Power-Path_Loss_dB + nakagami_fading_mean;
rx_power_alpha_mu_mean = Tx_Power-Path_Loss_dB + alpha_mu_fading_mean;



%% surf fading
range_quantization_par = 1;
range_iter = 0:range_quantization_par:400;

linear_fds_temp  = fds;
hist_range = min(linear_fds_temp(:,2)):max(linear_fds_temp(:,2));
linear_fds_temp(:,1) = floor(linear_fds_temp(:,1)/range_quantization_par)*range_quantization_par;
[range_mesh,bin_mesh] = meshgrid(range_iter,hist_range);
pdf_concat=[]
cmap_concat = []

for i=range_iter
    index = linear_fds_temp(:,1)==i;
    ddr = linear_fds_temp(index,:);
    a = size(ddr);

        [pdf,edge]=hist(ddr(:,2),hist_range);
        pdf = pdf/sum(pdf);
        cmap_temp  = pdf/max(pdf(:));
        pdf_concat = [pdf_concat,pdf'];
        cmap_concat = [cmap_concat,cmap_temp'];

%         title(['Range:' ,num2str(i)])
%         saveas(gcf,['../Matlab/Plots/',dataset,'Range:' ,num2str(i),'Fading Distribution.png'])
%         pause

end

h = figure
s=surf(range_mesh,bin_mesh,pdf_concat,cmap_concat,'FaceAlpha',1);
% s=surf(range_mesh,bin_mesh,pdf_concat,'FaceAlpha',1);
zlim([0,.3])
s.EdgeColor='none'
colorbar
title('Fading Process')
saveas(gcf,[input,'GT_Fading.fig'])
saveas(gcf,[input,'GT_Fading.png'])
%% nakagami Fading
range_quantization_par = 1;
range_iter = 0:range_quantization_par:400;

linear_fds_temp  = nakagami_dbm;
hist_range = min(linear_fds_temp(:,2)):max(linear_fds_temp(:,2));
linear_fds_temp(:,1) = floor(linear_fds_temp(:,1)/range_quantization_par)*range_quantization_par;
[range_mesh,bin_mesh] = meshgrid(range_iter,hist_range);
pdf_concat=[]
cmap_concat = []

for i=range_iter
    index = linear_fds_temp(:,1)==i;
    ddr = linear_fds_temp(index,:);
    a = size(ddr);

        [pdf,edge]=hist(ddr(:,2),hist_range);
        pdf = pdf/sum(pdf);
        cmap_temp  = pdf/max(pdf(:));
        pdf_concat = [pdf_concat,pdf'];
        cmap_concat = [cmap_concat,cmap_temp'];

%         title(['Range:' ,num2str(i)])
%         saveas(gcf,['../Matlab/Plots/',dataset,'Range:' ,num2str(i),'Fading Distribution.png'])
%         pause

end

figure
s=surf(range_mesh,bin_mesh,pdf_concat,cmap_concat,'FaceAlpha',1);
% s=surf(range_mesh,bin_mesh,pdf_concat,'FaceAlpha',1);
zlim([0,.3])
s.EdgeColor='none'
colorbar
title('Nakagami Fading Process')
saveas(gcf,[input,'nakagami_Fading.fig'])
saveas(gcf,[input,'nakagami_Fading.png'])
%% Alpha mu Fading
range_quantization_par = 1;
range_iter = 0:range_quantization_par:400;

linear_fds_temp  = alpha_mu_dbm;
hist_range = min(linear_fds_temp(:,2)):max(linear_fds_temp(:,2));
linear_fds_temp(:,1) = floor(linear_fds_temp(:,1)/range_quantization_par)*range_quantization_par;
[range_mesh,bin_mesh] = meshgrid(range_iter,hist_range);
pdf_concat=[]
cmap_concat = []

for i=range_iter
    index = linear_fds_temp(:,1)==i;
    ddr = linear_fds_temp(index,:);
    a = size(ddr);

        [pdf,edge]=hist(ddr(:,2),hist_range);
        pdf = pdf/sum(pdf);
        cmap_temp  = pdf/max(pdf(:));
        pdf_concat = [pdf_concat,pdf'];
        cmap_concat = [cmap_concat,cmap_temp'];

%         title(['Range:' ,num2str(i)])
%         saveas(gcf,['../Matlab/Plots/',dataset,'Range:' ,num2str(i),'Fading Distribution.png'])
%         pause

end

figure
s=surf(range_mesh,bin_mesh,pdf_concat,cmap_concat,'FaceAlpha',1);
% s=surf(range_mesh,bin_mesh,pdf_concat,'FaceAlpha',1);
zlim([0,.3])
s.EdgeColor='none'
colorbar
title('Alpha MU Fading Process')
saveas(gcf,[input,'Alpha_MU_Fading.fig'])
saveas(gcf,[input,'Alpha_MU_Fading.jpg'])
%% median plotting

% figure
% 
% lower = find(data_set(:,1)<=d_min-1,1,'last');
% upper = find(data_set(:,1)<=d_max,1,'last');
% ax1 = boxplot(data_set(:,2),data_set(:,1),'plotstyle','compact','colors','r','outliersize',10^-10);
% 
% title([input,'RSSI-Distance'])
% xlabel('Distance')
% ylabel('RSSI')
% ylim([-100 -40])
% saveas(gcf,[input,'Field-RSSI-Distance.png'])
% 
% figure
% ax2= boxplot(alpha_mu_recons(:,2),alpha_mu_recons(:,1),'plotstyle','compact','colors','g','outliersize',10^-10,'colors','g');
% title(['alpha mu RSSI-Distance'])
% xlabel('Distance')
% ylabel('RSSI')
% ylim([-100 -40])
% saveas(gcf,[input,'alpha-mu-RSSI-Distance.png'])
% 
% figure 
% boxplot(nakagami_recons(:,2),nakagami_recons(:,1),'plotstyle','compact','colors','b','outliersize',10^-10,'colors','b');
% title(['Nakagami RSSI-Distance'])
% xlabel('Distance')
% ylabel('RSSI')
% ylim([-100 -40])
% saveas(gcf,[input,'Nakagami-RSSI-Distance.png'])

%% fading box plot
am_trunc_index = alpha_mu_recons(:,2)>=-97;
am_trunc_dbm = alpha_mu_recons(am_trunc_index,:);
nak_trunc_index = nakagami_recons(:,2)>=-97;
nak_trunc_dbm = nakagami_recons(nak_trunc_index,:);
bin_size = 10 ;
scrsz = get(groot,'ScreenSize');
figure('Position',[1 1 1920 1080])

% ylim([-100 -40])

lower = find(data_set(:,1)<=d_min-1,1,'last');
upper = find(data_set(:,1)<=d_max,1,'last');
ax1 = boxplot(data_set(:,2),floor(data_set(:,1)/bin_size)*bin_size,'plotstyle','compact','colors','r','outliersize',10^-10);
hold on
ax2= boxplot(am_trunc_dbm(:,2),floor(am_trunc_dbm(:,1)/bin_size)*bin_size,'plotstyle','compact','colors','g','outliersize',10^-10,'colors','g');
xlabel('Distance')
ylabel('RSSI')
title([input,'\alpha-\mu(g) vs Field(r) RSSI-Distance'])
hold off
saveas(gcf,[input,'alpha-mu vs Field RSSI-Distance.png'])

figure('Position',[1 1 1920 1080])



ax1 = boxplot(data_set(:,2),floor(data_set(:,1)/bin_size)*bin_size,'plotstyle','compact','colors','r','outliersize',10^-10);
hold on
boxplot(nak_trunc_dbm(:,2),floor(nak_trunc_dbm(:,1)/bin_size)*bin_size,'plotstyle','compact','colors','b','outliersize',10^-10,'colors','b');
xlabel('Distance')
ylabel('RSSI')
ylim([-100 -40])
title([input,'Nakagami(b) vs Field(r) Fading-Distance'])
saveas(gcf,[input,'Nakagami vs Field RSSI-Distance.png'])
%% KLD PLOT
kld_alpha_mean = repmat(kld_alpha,size(kld_alpha,2),1);
kld_alpha_mean = tril(kld_alpha_mean,0);
kld_alpha_mean = sum(kld_alpha_mean,2)./(1:size(kld_alpha_mean,1))';
figure 
subplot(1,2,1)

plot(1:d_max,kld_alpha)
title('KLD \alpha-\mu')
subplot(1,2,2)
plot(1:d_max,kld_nak)
title('KLD Nakagami')
saveas(gcf,[input,'KLD.jpg'])
figure
subplot(1,2,1)

hval_cut = hval
hval_cut((hval_cut>1))=nan
hval_cut(hval_cut<1)=nan
plot(1:size(hval),hval_cut)
title('alpha-mu KS-Test Decision')
subplot(1,2,2)
plot(1:size(hval),pval)
title('alpha-mu KS-Test Asymptotic parameter')
saveas(gcf,[input,'alpha-mu RSSI-Distance.png'])

figure
subplot(1,2,1)
hval_cut = Nakhval
hval_cut(hval_cut>1)=nan
hval_cut(hval_cut<1)=nan

plot(1:size(hval_cut),hval_cut,'lineWidth',.1)
title('nakagami KS-Test Decision')
subplot(1,2,2)
plot(1:size(hval),Nakpval)
title('nakagami KS-Test Asymptotic parameter')
saveas(gcf,[input,'nakagami RSSI-Distance.png'])

figure
plot(1:size(rx_power_nakagami_mean),rx_power_nakagami_mean,1:size(rx_power_alpha_mu_mean),rx_power_alpha_mu_mean,1:size(rx_power_alpha_mu_mean),RSSI_per_meter_mean)
legend('nakagami','alpha mu','dataset')
title('Comparison of mean of each process')
saveas(gcf,[input,'Mean Comparison.png'])


