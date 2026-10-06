clc
clear

data=importdata('pmpulse.txt')
size_group=size(data);
N_group=size_group(1);

%% FUNCG: Define siganl
offset=0;%偏置
T_signal=1/(5E4); %s 序列长度
N_signal=10000; %points %采样数
N_pulse=70; %脉冲数
T_pulse=(3+0)*13.1728E-9; %脉冲间隔

dt_signal=T_signal/N_signal; %s  %每一个采样点长度
P_wave=round(T_pulse/dt_signal)                         %波形点数


for j=1:N_group
    
%% 自序列生成

a=data(j,:)
Body = zeros(P_wave*N_pulse,1);
for i=1:N_pulse
    pulse(i,1:P_wave)=a(i); 
end
for i=1:N_pulse
    Body((i-1)*P_wave+1:i*P_wave)=pulse(i,:);
end
Body=Body + offset;%加上偏置

%% 组成脉冲序列
space(1:N_signal-N_pulse*P_wave)=-.9;
waveform = [Body; space'];  

% % Plot the custom waveform to be generated

plot(waveform); 
axis([-0.1*N_signal 1.1*N_signal  -1.2 1.2])
set(gca, 'XTick', [0 N_signal])
set(gca, 'XTickLabel' ,{'0','10us'})
set(gca, 'YTick', [-1 -0.5 0 0.5 1])

b=num2str(j);
csvwrite(['C:\Users\HefeiLaptop\Desktop\tek_control\2022_6_10更新\tek_control - 副本\tek_control - 副本\tek_control - 副本\waveform_pm\waveform' b '.csv'],waveform);
pause(0.1);
end