%% 求取lora信号的csi相位

%% 1.读取数据
raw_data=readData('I:\Projects\files\LoRa_sense\phase_calculate\19430\TEST'); %读labview收的数据
% raw_data1=read_complex_binary('I:\Projects\files\LoRa_sense\phase_calculate\test01'); %读gnuradio数据
% 采样率
fs=1e6;
% 查看时域波形，找信号与噪声阈值
% plot(abs(raw_data1(1:fs/2)))

%% 提取信号部分
% (暂时没用到)
threshold=0.02; %设置信号与噪声区分的阈值
sign=raw_data(abs(raw_data)>threshold); 
% sign1=raw_data1(abs(raw_data1)>0.3);
% spectrogram(raw_data(1024*10.3:1024*22.8),64,60,64,fs,'yaxis','centered')

%% 2.产生标准chirp && preamble
Bw=125e3;
SF=7;
symbol_time = 2^SF/Bw;
t=1/fs:1/fs:symbol_time;
f1=-Bw/2+Bw/2/t(end)*t;
f2=Bw/2-Bw/2/t(end)*t;
% % % 自己产生
% % std_chirp=exp(1j*2*pi*f1.*t);   % upchirp
% % std_dwchp=exp(1j*2*pi*f2.*t);   % downchirp
% % cor_chirp=repmat(std_chirp,1,10);
% % cor_chirp=[cor_chirp std_dwchp std_dwchp std_dwchp(1:0.25*length(std_dwchp))];
% % % 自己产生version2
% % std_chirp2=cos(2*pi*f1.*t)+1j*cos(2*pi*f1.*t+pi/2);   % upchirp
% % std_dwchp2=cos(2*pi*f2.*t)+1j*sin(2*pi*f2.*t+pi/2);   % downchirp
% % cor_chirp2=repmat(std_chirp2,1,10);
% % cor_chirp2=[cor_chirp2 std_dwchp2 std_dwchp2 std_dwchp2(1:0.25*length(std_dwchp2))];
% % plot(abs(std_chirp2))
% 调用matlab chirp函数产生
% 标准chirp
chirpI=chirp(t,-Bw/2,symbol_time,Bw/2,'linear',90);
chirpQ=chirp(t,-Bw/2,symbol_time,Bw/2,'linear',0);
std_chirp=chirpI+1j*chirpQ;
std_dwchp=conj(std_chirp);
% 标准preamble
premb=repmat(std_chirp,1,10);
std_premb=[premb std_dwchp std_dwchp std_dwchp(1:0.25*length(std_dwchp))];   % 原始信号premble中
  % 第9、10个chirp是调制过的用以标识网络号的？？
% plot(abs(std_chirp))
% spectrogram(std_premb,64,60,64,fs,'yaxis','centered')

%% 3.匹配滤波/相关以截取整个preamble(粗同步)
% 整个preamble相关（未做归一化，归一化后是否更好？？）
[s_cor,lag]=xcorr(diff(phase(raw_data(1024*0+1:1024*300))),diff(phase(std_premb)));
% plot(abs(s_cor))
% 取峰值和对应的延迟
[h,index]=findpeaks(abs(s_cor),'MinPeakHeight',540); % 取大于阈值(540)的峰值和对应的index
lag_array=lag(index); % premble的延迟向量
% 截取整个preamble部分放入一个另外的矩阵
pream_matr=[];
for i=1:length(lag_array)
    pream_matr=[pream_matr;raw_data(lag_array(i)-3:lag_array(i)+1024*12+3)];   % 截取preamble放入另
        % 一矩阵,前后各多截取3个采样点，保证不会少截，3>=该段时间内采样时钟漂移
end
% spectrogram(pream_matr(1,:),64,60,64,'yaxis','centered')

%% 4.匹配滤波/相关以截取单个chirp(符号定时、细同步)
% 单个chirp符号定时与截取

% 做相关以求取符号定时
for i=1:size(pream_matr,1)
    [sig_chp_cor,sig_lag_arr]=xcorr(diff(phase(pream_matr(i,:))),diff(phase(std_chirp)));  % 与单个upchirp做相关
    [schp_cor_pk,index1]=findpeaks(sig_chp_cor,'MinPeakHeight',42);   % 找到峰值大于阈值42的点
    [sig_dchp_cor,sdchp_lag_arr]=xcorr(diff(phase(pream_matr(i,:))),diff(phase(std_dwchp)));  % 与单个downchirp做相关
    [sdchp_cor_pk,index2]=findpeaks(sig_dchp_cor,'MinPeakHeight',42);   % 找到峰值大于阈值42的点
    index_all=[index1 index2];
    for j=1:length(index_all)
        stru_sig_chp(i,j)={pream_matr(i,(sig_lag_arr(index_all(j)):sig_lag_arr(index_all(j))+1023))};  % 截取chirp放入另一cell中
    end
end

% scatter(1:length(sig_chp_cor),sig_chp_cor);hold on
% scatter(index1,schp_cor_pk);hold on
% scatter(1:length(sig_dchp_cor),sig_dchp_cor);hold on
% scatter(index2,sdchp_cor_pk);hold on
% spectrogram(cell2mat(stru_sig_chp(1,9)),64,60,64,fs,'yaxis','centered')