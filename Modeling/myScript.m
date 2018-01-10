clc
clearvars -except data_set
close all
clear
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

%% getting the mean

d_max = 899;
d_min = 4;

fixed_fading_distance = 1;

Tx_Power = 20;

index=[];

RSSI_per_meter_mean = zeros(d_max,1);
total_point_per_meter = zeros(d_max,1);
lower = find(data_set(:,1)<=d_min-1,1,'last');

for k=d_min:d_max
    
    upper = find(data_set(:,1)<=k,1,'last');
    RSSI_per_meter_mean(k,1)=median(data_set(lower+1:upper,2));
    total_point_per_meter(k) = upper - lower;
    lower = upper;
end

clear lower upper

% plot(index,RSSI_per_meter_mean(index),'.')

%%
param = struct('alpha',1,'eps',1,'tx_h',1,'freq',1)
distance= d_min:d_max;
alpha_list = 2.0:.1:2.01;
Epsilon_r_list = 1:.001:1.01;
CarrierFreq=5.89*10^9;
Tx_height = 1.61;
freq_list = linspace(CarrierFreq,CarrierFreq+10^9,10)
height_list = linspace(Tx_height,Tx_height+1,10)

Rx_height = Tx_height;

Speedof_Light=3*10^8;
alpha =2
eps = 1.01
CarrierFreq=5.89*10^9;
Tx_height = 1.61;
param.alpha = alpha;param.eps =eps;param.freq = CarrierFreq;param.tx_h=Tx_height
%% removing the path-loss
% for alpha = alpha_list
%     for Epsilon_r = Epsilon_r_list
%         
%         Path_Loss_dB = pl_estimator(param,d_min,d_max);
%         lower = find(data_set(:,1)<=d_min-1,1,'last');
% 
%         figure;
%         plot(d_min:d_max,Tx_Power-Path_Loss_dB(d_min:d_max),'r')
%         hold on 
%         plot(1:size(RSSI_per_meter_mean(:)),RSSI_per_meter_mean,'.b')
%         title(['alpha:',num2str(alpha),' e:',num2str(Epsilon_r)])
%     end
% end
%% Param Iterate
alpha =2
eps = 1.01
CarrierFreq=5.265*10^9;

Tx_height = 1.61;
Tx_height = 1.56
param.alpha = alpha;param.eps =eps;param.freq = CarrierFreq;param.tx_h=Tx_height

% i=1
% figure;
% for freq = linspace(CarrierFreq,CarrierFreq+10^8,9)
%     
%     param.freq= freq;
%     Path_Loss_dB = pl_estimator(param,d_min,d_max);
%     lower = find(data_set(:,1)<=d_min-1,1,'last');
% 
%     
%     subplot(3,3,i)
%     plot(d_min:d_max,Tx_Power-Path_Loss_dB(d_min:d_max),'r')
%     hold on 
%     plot(1:size(RSSI_per_meter_mean(:)),RSSI_per_meter_mean,'.b')
%     title(evalc(['disp(param)']))
%     i=i+1
% end
% param.alpha = alpha;param.eps =eps;param.freq = CarrierFreq;param.tx_h=Tx_height
% figure;
% i=1
% for t_h = linspace(Tx_height-.1,Tx_height+.1,9)
%     
%     param.tx_h= t_h;
%     Path_Loss_dB = pl_estimator(param,d_min,d_max);
%     lower = find(data_set(:,1)<=d_min-1,1,'last');
% 
%     
%     subplot(3,3,i);
%     plot(d_min:d_max,Tx_Power-Path_Loss_dB(d_min:d_max),'r');
%     hold on 
%     plot(1:size(RSSI_per_meter_mean(:)),RSSI_per_meter_mean,'.b')
%     title(evalc(['disp(param)']))
%     i=i+1
% end
% param.alpha = alpha;param.eps =eps;param.freq = CarrierFreq;param.tx_h=Tx_height
% figure;
% i=1
% for alpha_i = linspace(alpha,alpha+.1,9)
%     
%     param.alpha= alpha_i;
%     Path_Loss_dB = pl_estimator(param,d_min,d_max);
%     lower = find(data_set(:,1)<=d_min-1,1,'last');
% 
%     
%     subplot(3,3,i);
%     plot(d_min:d_max,Tx_Power-Path_Loss_dB(d_min:d_max),'r');
%     hold on 
%     plot(1:size(RSSI_per_meter_mean(:)),RSSI_per_meter_mean,'.b');
%     title(evalc(['disp(param)']));
%     i=i+1
% end
% param.alpha = alpha;param.eps =eps;param.freq = CarrierFreq;param.tx_h=Tx_height
% figure;
% i=1
% for eps = linspace(eps,eps+.1,9)
%     
%     param.eps= eps;
%     Path_Loss_dB = pl_estimator(param,d_min,d_max);
%     lower = find(data_set(:,1)<=d_min-1,1,'last');
% 
%     
%     subplot(3,3,i)
%     plot(d_min:d_max,Tx_Power-Path_Loss_dB(d_min:d_max),'r')
%     hold on 
%     plot(1:size(RSSI_per_meter_mean(:)),RSSI_per_meter_mean,'.b')
%     title(evalc(['disp(param)']))
%     i=i+1
% end
%% Grid search
loss_res = inf
alpha =2.0576
eps = 1.0036
% eps=1
CarrierFreq=5.89*10^9;

Tx_height = 1.61;
Tx_height = 1.4787
t_h = Tx_height
freq = CarrierFreq
%% params for freq 5.89 density 0-12
alpha=2.0108;
eps = 1.0043;
CarrierFreq=5.89*10^9;
Tx_height = 1.4787;

%% params for freq 5.89 density 12-130
% alpha=2.0237
% eps = 1.0091
% CarrierFreq=5.89*10^9;
% Tx_height = 1.4787

%% params for freq 5.89 density 30-50
% alpha=2.0576
% eps = 1.0036
% CarrierFreq=5.89*10^9;
% Tx_height = 1.4787
% for eps = linspace(eps,eps+.001,1000)
%     
% for alpha_i = linspace(alpha,alpha,1)
    
% for t_h = linspace(Tx_height-.5,Tx_height+.5,15000)
%     t_h
% for freq = linspace(CarrierFreq-10^5,CarrierFreq+10^5,20)
    th = Tx_height
    param.eps= eps;
    param.alpha= alpha;
    param.tx_h= t_h;
    param.freq= freq;
    Path_Loss_dB = pl_estimator(param,d_min,d_max);
    
    loss_temp = mean(abs(Tx_Power-Path_Loss_dB(1:100)-RSSI_per_meter_mean(1:100)).^2,'omitnan');
    if loss_temp<loss_res
        path_loss_res = Path_Loss_dB;
        param_res = param;
        loss_temp
        
        loss_res=loss_temp;
    end
    lower = find(data_set(:,1)<=d_min-1,1,'last');

% end
% end
% end
% end
figure
    plot(d_min:d_max,Tx_Power-path_loss_res(d_min:d_max),'r')
    hold on 
    plot(1:size(RSSI_per_meter_mean(:)),RSSI_per_meter_mean,'.b')
    title(evalc(['disp(param_res)']))
%% Chosen params
alpha=2.02
eps = 1.0148
CarrierFreq=5.2649*10^9;
Tx_height = 1.565

%% params for freq 5.89 density 0-12
alpha=2.0108;
eps = 1.0043;
CarrierFreq=5.89*10^9;
Tx_height = 1.4787;

%% params for freq 5.89 density 12-130
alpha=2.0237
eps = 1.0091
CarrierFreq=5.89*10^9;
Tx_height = 1.4787

%% params for freq 5.89 density 30-50
alpha=2.0576
eps = 1.0036
CarrierFreq=5.89*10^9;
Tx_height = 1.4787
