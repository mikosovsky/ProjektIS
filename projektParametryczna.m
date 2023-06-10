clear all;
clc;
clf
close all;
vec = load('robot_arm.dat');
Tp = 0.01;
fs = 1/Tp;
u = vec(:,1);
y = vec(:,2);
n = length(u);
M = n;
t=0:Tp:(M-1)*Tp;
percentProbekDoWyznaczenia = 0.7;
%% Wybrać rząd mianownika i licznika
mianownik=5;
licznik=1;
dlugoscMacierzy = n-max(licznik,mianownik);
Phi = zeros(dlugoscMacierzy,licznik+mianownik);
%% Pętla wypełniająca macierz Phi
for i=1:dlugoscMacierzy
    for j=1:mianownik
        Phi(i,j) = -y(i+max(mianownik,licznik)-j);
    end
    for j=1:licznik
        Phi(i,mianownik+j) = u(i+max(licznik,mianownik)-j+1);
    end
end
%% Obliczanie wektora theta
y_n = y(max(mianownik,licznik)+1:n);
y_n1 = y_n(1:length(y_n)*percentProbekDoWyznaczenia);
Phi1 = Phi(1:length(Phi)*percentProbekDoWyznaczenia,:);
theta = (Phi1'*Phi1)^(-1)*Phi1'*y_n1;
%% Wektor błedów predykcji i kwadratowy wskaźnik jakości
epsilonN = y_n - Phi*theta;
Vn = sum(epsilonN.^2)
Gz = tf(theta(licznik+mianownik:-1:mianownik+1)',theta(mianownik:-1:1)',0.01);
bode(Gz)
%% Wektor do weryfikacji ostatnich 30% zbioru
Phi2 = Phi(length(Phi)*percentProbekDoWyznaczenia:end,:);
y_test = Phi2*theta;
y_r = zeros(floor(n*(1-percentProbekDoWyznaczenia))-4,1);
y_g = y(floor(n*percentProbekDoWyznaczenia)+max(licznik,mianownik):end);
t_g = t(1:floor(n*(1-percentProbekDoWyznaczenia))-4);
y_r(1:end) = y_test(4:end);
y_g = y_g(1:length(y_r));
figure
plot(t_g,y_g,'b')
hold on
plot(t_g,y_r,'r--')
legend("Sygnał oryginalny","30% próbek nie branych pod uwage w MNK")
title("Porównanie oryginalnego sygnału z sygnałem estymowanym")
xlabel("nT_p [s]")
epsilonN = y_g - y_r;
Vn = sum(epsilonN.^2)

%% Porównanie rzeczywistych wartości z naszym obiektem
y_ans = y(length(y)*percentProbekDoWyznaczenia:end);
y_test = Phi*theta;
y_r = zeros(n,1);
y_r(max(licznik,mianownik)+1:end) = y_test(1:end);
figure
plot(t(1:end-max(mianownik,licznik)),y_test,'r--')
hold on
plot(t,y,'b')
legend("LS bez opóźnienia","Sygnał oryginalny")
title("Porównanie oryginalnego sygnału z sygnałem estymowanym bez opóźnienia")
xlabel("nT_p [s]")
figure
plot(t,y,'b')
hold on
plot(t,y_r,'r--')
legend("Sygnał oryginalny","LS z opóźnieniem")
title("Porównanie oryginalnego sygnału z sygnałem estymowanym z opóźnieniem")
xlabel("nT_p [s]")
