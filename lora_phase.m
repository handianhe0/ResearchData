%% ��ȡlora�źŵ�csi��λ

%% 1.��ȡ����
raw_data=readData('I:\Projects\files\LoRa_sense\phase_calculate\19430\TEST'); %��labview�յ�����
% raw_data1=read_complex_binary('I:\Projects\files\LoRa_sense\phase_calculate\test01'); %��gnuradio����
% ������
fs=1e6;
% �鿴ʱ���Σ����ź���������ֵ
% plot(abs(raw_data1(1:fs/2)))

%% ��ȡ�źŲ���
% (��ʱû�õ�)
threshold=0.02; %�����ź����������ֵ���ֵ
sign=raw_data(abs(raw_data)>threshold); 
% sign1=raw_data1(abs(raw_data1)>0.3);
% spectrogram(raw_data(1024*10.3:1024*22.8),64,60,64,fs,'yaxis','centered')

%% 2.������׼chirp && preamble
Bw=125e3;
SF=7;
symbol_time = 2^SF/Bw;
t=1/fs:1/fs:symbol_time;
f1=-Bw/2+Bw/2/t(end)*t;
f2=Bw/2-Bw/2/t(end)*t;
% % % �Լ�����
% % std_chirp=exp(1j*2*pi*f1.*t);   % upchirp
% % std_dwchp=exp(1j*2*pi*f2.*t);   % downchirp
% % cor_chirp=repmat(std_chirp,1,10);
% % cor_chirp=[cor_chirp std_dwchp std_dwchp std_dwchp(1:0.25*length(std_dwchp))];
% % % �Լ�����version2
% % std_chirp2=cos(2*pi*f1.*t)+1j*cos(2*pi*f1.*t+pi/2);   % upchirp
% % std_dwchp2=cos(2*pi*f2.*t)+1j*sin(2*pi*f2.*t+pi/2);   % downchirp
% % cor_chirp2=repmat(std_chirp2,1,10);
% % cor_chirp2=[cor_chirp2 std_dwchp2 std_dwchp2 std_dwchp2(1:0.25*length(std_dwchp2))];
% % plot(abs(std_chirp2))
% ����matlab chirp��������
% ��׼chirp
chirpI=chirp(t,-Bw/2,symbol_time,Bw/2,'linear',90);
chirpQ=chirp(t,-Bw/2,symbol_time,Bw/2,'linear',0);
std_chirp=chirpI+1j*chirpQ;
std_dwchp=conj(std_chirp);
% ��׼preamble
premb=repmat(std_chirp,1,10);
std_premb=[premb std_dwchp std_dwchp std_dwchp(1:0.25*length(std_dwchp))];   % ԭʼ�ź�premble��
  % ��9��10��chirp�ǵ��ƹ������Ա�ʶ����ŵģ���
% plot(abs(std_chirp))
% spectrogram(std_premb,64,60,64,fs,'yaxis','centered')

%% 3.ƥ���˲�/����Խ�ȡ����preamble(��ͬ��)
% ����preamble��أ�δ����һ������һ�����Ƿ���ã�����
[s_cor,lag]=xcorr(diff(phase(raw_data(1024*0+1:1024*300))),diff(phase(std_premb)));
% plot(abs(s_cor))
% ȡ��ֵ�Ͷ�Ӧ���ӳ�
[h,index]=findpeaks(abs(s_cor),'MinPeakHeight',540); % ȡ������ֵ(540)�ķ�ֵ�Ͷ�Ӧ��index
lag_array=lag(index); % premble���ӳ�����
% ��ȡ����preamble���ַ���һ������ľ���
pream_matr=[];
for i=1:length(lag_array)
    pream_matr=[pream_matr;raw_data(lag_array(i)-3:lag_array(i)+1024*12+3)];   % ��ȡpreamble������
        % һ����,ǰ������ȡ3�������㣬��֤�����ٽأ�3>=�ö�ʱ���ڲ���ʱ��Ư��
end
% spectrogram(pream_matr(1,:),64,60,64,'yaxis','centered')

%% 4.ƥ���˲�/����Խ�ȡ����chirp(���Ŷ�ʱ��ϸͬ��)
% ����chirp���Ŷ�ʱ���ȡ

% ���������ȡ���Ŷ�ʱ
for i=1:size(pream_matr,1)
    [sig_chp_cor,sig_lag_arr]=xcorr(diff(phase(pream_matr(i,:))),diff(phase(std_chirp)));  % �뵥��upchirp�����
    [schp_cor_pk,index1]=findpeaks(sig_chp_cor,'MinPeakHeight',42);   % �ҵ���ֵ������ֵ42�ĵ�
    [sig_dchp_cor,sdchp_lag_arr]=xcorr(diff(phase(pream_matr(i,:))),diff(phase(std_dwchp)));  % �뵥��downchirp�����
    [sdchp_cor_pk,index2]=findpeaks(sig_dchp_cor,'MinPeakHeight',42);   % �ҵ���ֵ������ֵ42�ĵ�
    index_all=[index1 index2];
    for j=1:length(index_all)
        stru_sig_chp(i,j)={pream_matr(i,(sig_lag_arr(index_all(j)):sig_lag_arr(index_all(j))+1023))};  % ��ȡchirp������һcell��
    end
end

% scatter(1:length(sig_chp_cor),sig_chp_cor);hold on
% scatter(index1,schp_cor_pk);hold on
% scatter(1:length(sig_dchp_cor),sig_dchp_cor);hold on
% scatter(index2,sdchp_cor_pk);hold on
% spectrogram(cell2mat(stru_sig_chp(1,9)),64,60,64,fs,'yaxis','centered')