clear all;
clc;
clf
close all;
vec = load('robot_arm.dat');
% vec = load('cstr.dat');
Tp = 0.01;
fs = 1/Tp;
u = vec(:,1);
y = vec(:,2);
n = length(u);
M = n;
t=0:Tp:(M-1)*Tp;

%% SSYGNAŁ NA WEJŚCIU I WYJŚCIU
figure
plot(u)
title("Sygnał wejściowy")
figure
plot(y)
title("Sygnał wyjściowy")


%% KORELACJE OBCIĄŻONE
ruy = xcorr(y,u,'biased');
ruu = xcorr(u,u,'biased');
ryy = xcorr(y,y,'biased');
figure
xcorrPlotLen = floor(length(ruy)/2);
plot(-xcorrPlotLen:xcorrPlotLen,ruy)
title("Wykres korelacji wejścia i wyjścia")
xlim([-xcorrPlotLen,xcorrPlotLen])
figure
plot(-xcorrPlotLen:xcorrPlotLen,ruu)
title("Wykres autokorelacji wejścia")
xlim([-xcorrPlotLen,xcorrPlotLen])
figure
plot(-xcorrPlotLen:xcorrPlotLen, ryy)
title("Wykres autokorelacji wyjścia")
xlim([-xcorrPlotLen,xcorrPlotLen])


%% ANALIZA KORELACJI I SPRAWDZANIE MACIERZY Ruu
Ruu = zeros(M);
ruy = ruy(n:end);
ruu = ruu(n:end);
for i=1:M
    Ruu(i,i:M) = ruu(1:M-i+1);
    Ruu(i:M,i) = ruu(1:M-i+1);
end
det(Ruu)
gM = -pinv(Ruu)*ruy(1:M); % pseudoodwrotna macierz
figure
plot(t,gM)
title("Odpowiedź impulsowa")
h = zeros(1,M);
for i=1:M
    h(i) = sum(gM(1:i));
end
figure
plot(t,h)
title("Odpowiedź skokowa")
% Te wykresy nie mają sensu ponieważ det(Ruu) = 0, więc jak można zauważyć
% odpowiedź skokowa i impulsowa oscylują

%% GĘSTOŚĆ WIDMOWA MOCY SYGNAŁU WEJŚCIOWEGO
[Puu, f] = periodogram(u,[],n,fs);
figure
plot(f,10*log10(Puu))
title("Wykres analizy widmowej mocy")
ylabel("Widmo mocy [dB]")
xlabel("Częstotliwość [Hz]")
[Pyy, f] = periodogram(y,[],n,fs);
figure
plot(f,10*log10(Pyy))
title("Wykres analizy widmowej mocy Pyy")
ylabel("Widmo mocy [dB]")
xlabel("Częstotliwość [Hz]")
figure
plot(f,20*log10(Pyy./Puu))
title("Wykres analizy widmowej mocy Pyy/Puu")
ylabel("Widmo mocy [dB]")
xlabel("Częstotliwość [Hz]")
nfft = 256;
window = hanning(1024);
noverlap = 128;
dflag = 'none';
figure
[Puu, f] = periodogram(u,window,nfft,fs);
plot(f,10*log10(Puu))
title("Wykres analizy widmowej mocy z oknem Hanning'a")
ylabel("Widmo mocy [dB]")
xlabel("Częstotliwość [Hz]")
[Pyy, f] = periodogram(y,window,nfft,fs);
figure
plot(f,10*log10(Puu))
title("Wykres analizy widmowej mocy z oknem Hanning'a Pyy")
ylabel("Widmo mocy [dB]")
xlabel("Częstotliwość [Hz]")
%% ZWIĄZEK PARSEVALA
E_time = sum(u.^2)/n
E_freq = sum(Puu)/length(Puu)


%% WZAJEMNA GĘSTOŚĆ
Puy = csd(u,y,nfft,window',noverlap,dflag);
figure
plot(f,10*log10(Puy))
ylabel("Widmo mocy [dB]")
xlabel("Częstotliwość [Hz]")
title("Wykres wzajemnej gęstości widmowej")


%% PRZEDZIAŁY UFNOŚCI
p = 0.5;
k = fix((n-noverlap)/(length(window)-noverlap));
alfa = 1-p;
confid = 2*k*Puu*(1./chi2inv([1-alfa/2 alfa/2],2*k));

%% ESTYMATOR FUNKCJI PRZEJŚCIA
h = ones(1,10)/10;
yTest = filter(h,1,u);
[Hest, f] = tfe(u,y,256,fs,256,128,'none');
H=freqz(h,1,f,fs);
figure
plot(f,abs(H))
title("Rzeczywista amplituda funkcji przejścia")
figure
plot(f,abs(Hest));
title("Estymata amplitudy funkcji przejścia")