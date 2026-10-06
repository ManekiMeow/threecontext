clear
clc
Time_rand=round(rand*10000);
%% Connect to the afg3022C using visa-usb connection 2019.11.19
%  Having TekVisa installed.
%  For finding the visa address, check matlab-APPS-InstrumentControlBox-VISA-USB
%  To know how to connect the afg, always check 'Session Log'
%  Inside the InstrumentControlBox-VISA-USB
VisaObjFunCG = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0356::C013812::0::INSTR', 'Tag', '');
%USB串口名输入及串口对象建立
if isempty(VisaObjFunCG)
    VisaObjFunCG = visa('KEYSIGHT', 'USB0::0x0699::0x0356::C013812::0::INSTR');
else
    fclose(VisaObjFunCG);
    VisaObjFunCG = VisaObjFunCG(1);
end
% Set buffersize before generate the waveform whlie afg connected but not open
% Cannot set buffersize when afg is open
VisaObjFunCG.OutputBufferSize = 10240000;
% Open function generator
fopen(VisaObjFunCG);
%% FUNCG: Reset the function generator
 fprintf(VisaObjFunCG, '*RST');
 fprintf(VisaObjFunCG, '*CLS'); 
 
%% Connect to MDO3024
N_time=10;
Lenth_pulsegroup=6500;

visa_brand = 'ni';
visa_address = 'USB0::0x0699::0x0408::C060964::0::INSTR';  %USB串口号输入
%'USB0::0x0699::0x0408::C060964::0::INSTR'
buffer = 60 * 1024; %20 KiB
record = 10000;
waveout=zeros(record,1);

DpoMsoMdo = instrfind('Type', 'visa-tcpip', 'RsrcName', visa_address , 'Tag', '');
%创建串口通信对象
if isempty(DpoMsoMdo)
    DpoMsoMdo = visa(visa_brand, visa_address, 'InputBuffer', buffer, ...
    'OutputBuffer', buffer);
else
    fclose(DpoMsoMdo);
    DpoMsoMdo = DpoMsoMdo(1);
end

% Connect to instrument object
fopen(DpoMsoMdo);
 
 
%% FUNCG: Define siganl

N_group=40;
delay=0.1;%延迟
tic
for i=278:277+N_group
    b=num2str(i);
   waveform=csvread(['C:\Users\HefeiLaptop\Desktop\tek_control\2022_6_10更新\tek_control - 副本\tek_control - 副本\tek_control - 副本\waveform_im\waveform' b '.csv']);
   waveform_pm=csvread(['C:\Users\HefeiLaptop\Desktop\tek_control\2022_6_10更新\tek_control - 副本\tek_control - 副本\tek_control - 副本\waveform_pm\waveform' b '.csv']); 
%% 画出图像
N_size=size(waveform);    
N_signal=N_size(1);  
waveform=waveform;
subplot(221)
plot(waveform); 
axis([-0.1*N_signal 1.1*N_signal  -0.5 1.2])
set(gca, 'XTick', [0 N_signal 2*N_signal])
set(gca, 'XTickLabel' ,{'0','2000ns','4000ns'})
set(gca, 'YTick', [0 0.5 1])


% waveform_pm=waveform_pm;
subplot(222)
plot(waveform_pm); 
axis([-0.1*N_signal 1.1*N_signal  -0.5 1.2])
set(gca, 'XTick', [0 N_signal 2*N_signal])
set(gca, 'XTickLabel' ,{'0','2000ns','4000ns'})
set(gca, 'YTick', [0 0.5 1])
% Normalize waveform
%waveform = waveform ./ (max(waveform)); 


%% FUNCG: Convert waveform 
%  AFG3022C is a 14 bits device
%  Convert the double values integer values between 0 and 16382
%  As required by the instrument
waveform =  round((waveform + 1.0)*8191);
waveformLength = length(waveform);
%  Encode variable 'waveform' into binary waveform data for AFG.
%  programmer manual for bit definitions.
binblock = zeros(2 * waveformLength, 1);
binblock(2:2:end) = bitand(waveform, 255);
binblock(1:2:end) = bitshift(waveform, -8);
binblock = binblock';
% Build binary block header
bytes = num2str(length(binblock));
header = ['#' num2str(length(bytes)) bytes];

%% waveform_pm
waveform_pm =  round((waveform_pm + 1.0)*8191);
waveformLength = length(waveform_pm);
%  Encode variable 'waveform' into binary waveform data for AFG.
%  programmer manual for bit definitions.
binblock_pm = zeros(2 * waveformLength, 1);
binblock_pm(2:2:end) = bitand(waveform_pm, 255);
binblock_pm(1:2:end) = bitshift(waveform_pm, -8);
binblock_pm = binblock_pm';
% Build binary block header
bytes_pm = num2str(length(binblock_pm));
header_pm = ['#' num2str(length(bytes_pm)) bytes_pm];


%% FUNCG: Major settings for Function generator (tek afg3022C)
% Resets the contents of edit memory and define the length of signal
% Transfer the custom waveform from MATLAB to edit memory of instrument
% Set the source to EXTernal otherwise the generator will keep triggering
% Force a trigger event immediately
% For specific settings, please check the AFG3000 series programmer manual

fprintf(VisaObjFunCG, ['DATA:DEF EMEM1, ' N_signal ';']); %1001
fwrite(VisaObjFunCG, [':TRACE EMEM1, ' header binblock ';'], 'uint8');
fprintf(VisaObjFunCG, 'SOURCE1:BURSt:STATe ON');          %1通道建立
fprintf(VisaObjFunCG, 'SOURCE1:BURSt:MODE TRIGgered');     %波形模式
%fprintf(VisaObjFunCG, 'TRIGger:SEQuence:SOURce EXTernal');
fprintf(VisaObjFunCG, 'SOURCE1:BURSt:NCYCles 1');       %Burst重复的周期
fprintf(VisaObjFunCG, 'SOURCE1:BURSt:TDElay 200ns');       %trigger delay 的时间
fprintf(VisaObjFunCG, 'SOURCE1:FUNCTION EMEM'); 
fprintf(VisaObjFunCG, 'SOURCE1:FREQUENCY 5E4');                %1通道序列重复频率<1E6
fprintf(VisaObjFunCG, 'SOURCE1:VOLTAGE:AMPLITUDE 8.00');         %1通道振幅
fprintf(VisaObjFunCG, 'SOURCE1:VOLTAGE:OFFSET 0.70');          %1通道偏置
fprintf(VisaObjFunCG, ':OUTP1 ON');                            %1通道输出
fprintf(VisaObjFunCG, 'TRIGger:SEQuence:TIMer 100us');
pause(1)

fprintf(VisaObjFunCG, ['DATA:DEF EMEM2, ' N_signal ';']); %1001
fwrite(VisaObjFunCG, [':TRACE EMEM2, ' header_pm binblock_pm ';'], 'uint8');
fprintf(VisaObjFunCG, 'SOURCE2:BURSt:STATe ON');          %2通道建立
fprintf(VisaObjFunCG, 'SOURCE2:BURSt:MODE TRIGgered');
%fprintf(VisaObjFunCG, 'TRIGger:SEQuence:SOURce EXTernal');
fprintf(VisaObjFunCG, 'SOURCE2:BURSt:NCYCles 1');
fprintf(VisaObjFunCG, 'SOURCE2:BURSt:TDElay 210ns');       %trigger delay 的时间
fprintf(VisaObjFunCG, 'SOURCE2:FUNCTION EMEM2'); 
fprintf(VisaObjFunCG, 'SOURCE2:FREQUENCY 5E4');                %2通道序列重复频率<1E6
fprintf(VisaObjFunCG, 'SOURCE2:VOLTAGE:AMPLITUDE 9.32');         %2通道振幅
fprintf(VisaObjFunCG, 'SOURCE2:VOLTAGE:OFFSET 0 ');          %2通道偏置
fprintf(VisaObjFunCG, ':OUTP2 ON');                            %2通道输出
fprintf(VisaObjFunCG, 'TRIGger:SEQuence:TIMer 100us');

pause(delay)
%% read
for j=1:N_time
    
fwrite(DpoMsoMdo, 'curve?');  % main command
samples = fread(DpoMsoMdo, record, 'int8');
waveout=samples;

%% 导出数据
b=num2str(i);
C=num2str(j);
RR=num2str(Time_rand);
csvwrite(['data5\data' RR '___' b '_' C '.csv'],waveout);

end

%% 显示波形
subplot(223)
% plot samples
plot(waveout);
i
end
toc
%%
fclose(VisaObjFunCG);                                          %关闭串口，不关闭会影响下次打开
clear VisaObjFunCG;
fclose(DpoMsoMdo); 
delete(DpoMsoMdo); 
clear DpoMsoMdo;
Time_rand