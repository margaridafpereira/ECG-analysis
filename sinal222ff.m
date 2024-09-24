clear all
close all
clc

%% Sinal

load("222.mat");

G = 200; % ganho(200 adu/mV)

ecg = val/G;

fc = input("Frequência de corte (50/60Hz) \n", 'i')
display("\n")
Fs = input("Frequência de amostragem \n")
display("\n")

% Carateristica do sinal
N = size(val,2);
t_total = (N/Fs);

% Representação em função do tempo
Ts = 1/Fs; % tempo de amostragem
t = 0 : Ts : (t_total-Ts);
##t_ms = t .* 1000;
subplot(2,4,1)
plot(t, val);
title('Sinal, s(n), corrompido com ruído');
xlabel('Tempo /s');
ylabel('Amplitude /mV');

% Representação em função da frequencia
NFFT=2^nextpow2(N); % próxima potência de 2 do comprimento de s
f = 0 : Fs/NFFT : Fs - Fs/NFFT;
S = fft(ecg, NFFT)*2/N;
subplot(2,4,5)
plot(f, abs(S));
axis([0 Fs/2]);
title('Amplitude Espectral de S em função da Frequência');
ylabel('|S| /mV');
xlabel('Frequência (Hz)');

% FIR passa alto 1 Hz
pkg load signal
ordem = 1000;
a = 1;
b = fir1(ordem,1/(Fs/2),'high');

% filtrar sinal
ecg_f = filtfilt(b,a,ecg);
subplot(2,4,2)
plot(t, ecg_f);
title('Filtro FIR Passa-alto 1 Hz');
xlabel('Tempo /s');
ylabel('Amplitude /mV');

S_f = fft(ecg_f, NFFT)*2/N;

subplot(2,4,6)
plot(f, abs(S_f));
axis([0 Fs/2]);
title('Amplitude Espectral de S em função da Frequência');
ylabel('|S| /mV');
xlabel('Frequência (Hz)');

% FIR filter to remove 60Hz from an ECG
pkg load signal
ordem = 200;
if fc = 60
  limi = 59;
  lims = 61;
elseif fc = 50
  limi = 49;
  lims = 51;
end

a = 1;
b = fir1(ordem,[limi lims]/(Fs/2),'stop');

% filtrar sinal
ecg_fi = filtfilt(b,a,ecg_f);

subplot(2,4,3)
plot(t, ecg_fi);
title('Filtro FIR Passa-banda 60 Hz');
xlabel('Tempo /s');
ylabel('Amplitude /mV');

S_f = fft(ecg_fi, NFFT)*2/N;

subplot(2,4,7)
plot(f, abs(S_f));
axis([0 Fs/2]);
title('Amplitude Espectral de S em função da Frequência');
ylabel('|S| /mV');
xlabel('Frequência (Hz)');


% Filtro Cheby2 passa baixo 65 Hz
Fc1 = 65;
n = 4;
Rs = 20; 

H = FiltroCheby2(Fs, Fc1, NFFT, Rs, n, 'baixo');
S_filt = H.*S_f;

s = ifft(S_filt);
s = s(1:N)*(N/2);

subplot(2,4,4);
plot(t, s);
title('Sinal filtrado pelo Filtro Chevy2 passa baixo');
xlabel('Tempo(s)');
ylabel('|S_filt| /mV');

subplot(2,4,8)
plot(f, abs(S_filt));
axis([0 Fs/2]);
title('Amplitude Espectral de S em função da Frequência');
ylabel('|S| /mV');
xlabel('Frequência (Hz)');

% Encontrar picos onda R

ecg = s;

figure
hold on;
plot(t,ecg);
xlabel('Tempo /s');
ylabel('Amplitude /mV');

tmin = 0.5; % minima distância RR é 0.6
Nmin = tmin*Fs;

[pks, idx] = findpeaks(ecg,"DoubleSided","MinPeakHeight",0.4,"MinPeakDistance", Nmin);

scatter (t(idx), ecg(idx), "r");

% Encontrar picos onda P

amp_maior = 0;
ecg = real(ecg);
for i = 1 : length(idx) 
  % QRS tem duração entre 0,06 e 0,10 segundos
  % Intervalo PR tem duração entre entre 0,12 e 0,20 segundos
  % Distância entre picos máxima = 0,05+0,20 = 0,25s
  % Distância entre picos minima = 0,03+0,06 = 0,09s
  
  imin = 0.12*Fs; % primeiro indice maior que o intervalo PR
  imin = floor(imin); % arredonda para o inteiro mais próximo <= a esse elemento
  t1 = 0.25;
  t2 = 0.09;
  N1 = t1*Fs;
  N2 = t2*Fs;
  i_min = idx(i)-N1;
  i_min = floor(i_min); % arredonda para o inteiro mais próximo <= a esse elemento
  i_max = idx(i)-N2;
  i_min = ceil(i_min); % arredonda para o inteiro mais próximo >= a esse elemento
  max = -100; 

  if idx(1) > N1 
    for j = i_min : i_max
      if ecg(j) > max
        max = ecg(j);
        max_i = j;
      end
    end
    pks_p(i) = max;
    idx_p(i) = max_i;
  elseif idx(1) > imin
    if i == 1
      for j = 1 : i_max
        if ecg(j) > max
          max = ecg(j);
          max_i = j;
        end
      end
      pks_p(i) = max;
      idx_p(i) = max_i;        
    else
      for j = i_min : i_max
        if ecg(j) > max
          max = ecg(j);
          max_i = j;
        end
      end
      pks_p(i) = max;
      idx_p(i) = max_i;    
    end  
  else
    if i ~= 1
      for j = i_min : i_max
        if ecg(j) > max
          max = ecg(j);
          max_i = j;
        end
      end
      pks_p(i-1) = max;
      idx_p(i-1) = abs(max_i);    
    end
  end
  
end

display("\n")
scatter (t(idx_p), pks_p, "b");

% Duração ondas P

for i = 1 : length(idx_p)
  % duração média 0,11 segundos, arritmia poderá ser maior  
  tmin = 0.05;
  Nmin = t1*Fs;
  t1 = 0.08;
  N1 = t1*Fs;
  i_min = idx_p(i)-N1;
  i_min = floor(i_min); % arredonda para o inteiro mais próximo <= a esse elemento
  i_max = idx_p(i)+N1; 
  i_min = ceil(i_min); % arredonda para o inteiro mais próximo >= a esse elemento
  
  % inicio de cada onda P (supondo que é o valor minimo nesse intervalo)  
  if idx_p(i) > Nmin
    if idx_p(i) > N1     
      min = 10;
      min_i = 0;
      for j = i_min : idx_p(i)
        if ecg(j) < min
          min = ecg(j);
          min_i = j;
        end
      end
      ind_min(i) = min_i;
      Amin(i) = min;
    else    
      min = 10;
      min_i = 0;
      for j = 1 : idx_p(i)
        if ecg(j) < min
          min = ecg(j);
          min_i = j;
         end
      end
      ind_min(i) = min_i;
      Amin(i) = min;
    end
  else
    ind_min(i) = 1;
    Amin(i) = ecg(1);
  end

  % fim de cada onda P (supondo que é o valor minimo nesse intervalo)
  if (N-idx_p(i)) > Nmin
    if (N-idx_p(i)) > N1      
      min = 10;
      min_i = 0;
      for k = (idx_p(i)): i_max
        if ecg(k) < min
          min = ecg(k);
          min_i = k;
        end
      end
      ind_max(i) = min_i;
      Amax(i) = min;  
    else  
      min = 10;
      min_i = 0;
      for k = (idx_p(i)): N
        if ecg(k) < min
          min = ecg(k);
          min_i = k;
        end
      end
      ind_max(i) = min_i;
      Amax(i) = min;
    end
  else
    ind_max(i) = N;
    Amax(i) = ecg(N);
  end
       
  % duração ondas P
  t_max = t(ind_max(i));
  t_min = t(ind_min(i));
  dur_p(i) = t_max - t_min;
 
end

scatter (t(ind_min), Amin, "g");
scatter (t(ind_max), Amax, "g");

% Frequência cardíaca
n_ciclos = length(idx_p);
t_min = t_total/60;
freq_cardiaca = n_ciclos/t_min;

% Menu estudo de patologias associadas aos valores de duração e amplitude
% ondas P e frequência cardíaca

choice = menu('Estudo de patologia', 'Crescimento átrio direito', 'Crescimento átrio esquerdo','Ritmo cardiaco','Todas as anteriores')
switch choice
  case 1
    amp_maior = 0;
    for i = 1 : length(pks_p)
      if pks_p(i) > 0.25
        disp('Amplitude acima do normal no pico: '), disp (i)
        disp('Valor amplitude: '), disp (pks_p(i))
        amp_maior += 1;
      end  
    end 
    if amp_maior > 0
      disp('Possibilidade de existência de crescimento do átrio direito!!')
      disp('Nota: analisar amplitudes ondas P.')
      display("\n") 
    else
      disp('Não foi detetada nenhuma irregularidade.')
    end
  case 2
    dur_maior = 0;
    for i = 1 : length(dur_p)
      if dur_p(i) > 0.20
        disp('Duração acima do normal na onda: '), disp (i)
        disp('Valor: '), disp (dur_p(i))
        display("\n")
        dur_maior += 1;
      end
    end 
    if dur_maior > 0
      disp('Possibilidade de existência de crescimento do átrio esquerdo!')
      disp('Nota: analisar duração ondas P.')
      display("\n")
    else
      disp('Não foi detetada nenhuma irregularidade.')
    end 
  case 3
    disp('Frequência cardiaca: '), disp (freq_cardiaca)
    if freq_cardiaca < 60
        disp('Possível situação de bradicardia!')
        display("\n")
    elseif freq_cardiaca > 60 && freq_cardiaca < 100
        disp('Não foi detetada nenhuma irregularidade.')
        display("\n")
    elseif freq_cardiaca > 100 && freq_cardiaca < 150
        disp('Possível situação de taquicardia!')
        display("\n")
    elseif freq_cardiaca > 150 && freq_cardiaca < 200
        disp('Possível situação de taquicardia paroxística!')
        display("\n")
    elseif freq_cardiaca > 200 && freq_cardiaca < 250
        disp('Possível situação de taquicardia paroxística e/ou flutter!')
        display("\n")
    elseif freq_cardiaca > 250 && freq_cardiaca < 300
        disp('Possível situação de taquicardia!')
        display("\n")
    end 
 case 4
    amp_maior = 0;
    for i = 1 : length(pks_p)
      if pks_p(i) > 0.25
        disp('Amplitude acima do normal no pico: '), disp (i)
        disp('Valor amplitude: '), disp (pks_p(i))
        amp_maior += 1;
      end  
    end 
    
    dur_maior = 0;
    for i = 1 : length(dur_p)
      if dur_p(i) > 0.20
        disp('Duração acima do normal na onda: '), disp (i)
        disp('Valor: '), disp (dur_p(i))
        display("\n")
        dur_maior += 1;
      end
    end 
    
    if amp_maior > 0
      disp('Possibilidade de existência de crescimento do átrio direito!!')
      disp('Nota: analisar amplitudes ondas P.')
      display("\n")
    end
    if dur_maior > 0
      disp('Possibilidade de existência de crescimento do átrio esquerdo!')
      disp('Nota: analisar duração ondas P.')
      display("\n")
    end
    
    disp('Frequência cardiaca: '), disp (freq_cardiaca)
    if freq_cardiaca < 60
        disp('Possível situação de bradicardia!')
        display("\n")
    elseif freq_cardiaca > 100 && freq_cardiaca < 150
        disp('Possível situação de taquicardia!')
        display("\n")
    elseif freq_cardiaca > 150 && freq_cardiaca < 200
        disp('Possível situação de taquicardia paroxística!')
        display("\n")
    elseif freq_cardiaca > 200 && freq_cardiaca < 250
        disp('Possível situação de taquicardia paroxística e/ou flutter!')
        display("\n")
    elseif freq_cardiaca > 250 && freq_cardiaca < 300
        disp('Possível situação de taquicardia!')
        display("\n")
    end 
    
    if (freq_cardiaca > 60 && freq_cardiaca < 100) && dur_maior == 0 && amp_maior == 0
        disp('Não foi detetada nenhuma irregularidade.')
        display("\n")
    end
   
endswitch
  