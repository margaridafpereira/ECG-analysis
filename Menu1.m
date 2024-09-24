%amplitude espectral 1

clear all
close all
clc
%-------------------------------------------------------
pkg load signal

%% Sinal

load("222.mat");
%------------------------------------------------------
G = 200; % ganho(200 adu/mV)

ecg = val/G;

fc = input("Frequência de corte (50/60Hz) \n", 'i')
display("\n")
Fs = input("Frequência de amostragem \n")
display("\n")
Fc = 1; 

% Carateristica do sinal
N = size(val,2);
t_total = (N/Fs);


Ts = 1/Fs; % tempo de amostragem
t = 0 : Ts : (t_total-Ts);
t_ms = t .* 1000;


NFFT=2^nextpow2(N); % próxima potência de 2 do comprimento de s
f = 0 : Fs/NFFT : Fs - Fs/NFFT;
S = fft(ecg, NFFT)*2/N;

w = 0 : 2*pi/NFFT: 2*pi - 2*pi/NFFT;

%------------------------------------------------------------------------------------
##############################
#####Filtros digitais IIR#####
##############################
%------------------------------------------------------------------------
n = 4;#ordem do filtro 

choice = menu('Ação pretendida', 'Retirar componente dos 60 Hz', 'Eliminar ruído até aos 1 Hz','Remover o ruído provocado por frequências altas')
switch choice
  case 1
    %Sinal Original
    subplot(2,3,1);
    plot(t, val);
    title('Sinal, s(n), corrompido com ruído');
    xlabel('Tempo /s');
    ylabel('Amplitude /mV');
      
      
    %Filtro Butter
    H = FiltroButter(Fs, fc, NFFT, n, 'stopp', 0.25, 0.85);
    S_filt = H.*S;

    s = ifft(S_filt);
    s = s(1:N)*(N/2);

    subplot(2,3,2);
    plot(t, s);
    title('Sinal filtrado pelo Filtro Butter rejeita banda');
    xlabel('Tempo (s)');
    ylabel('|S_filt| /mV');
          
##    S_f = fft(s, NFFT)*2/N;
##
##    plot(f, abs(S_f));
##    axis([0 Fs/2]);
##    title('Amplitude Espectral de S em função da Frequência');
##    ylabel('|S| /mV');
##    xlabel('Frequência (Hz)');
    
    
    %Filtro Cheby1
    Rp = 5; 
    H1 = FiltroCheby1(Fs, fc, NFFT, Rp, n, 'stopp', 0.15, 0.85);
    S_filt1 = H1.*S;
    s1 = ifft(S_filt1);
    s1 = s1(1:N)*(N/2);

    subplot(2,3,3);
    plot(t, s1);
    title('Sinal filtrado pelo Filtro Cheby1 rejeita banda');
    xlabel('Tempo (s)');
    ylabel('|S_filt| /mV');
    
##    S_f1 = fft(s1, NFFT)*2/N;
##
##    plot(f, abs(S_f1));
##    axis([0 Fs/2]);
##    title('Amplitude Espectral de S em função da Frequência');
##    ylabel('|S| /mV');
##    xlabel('Frequência (Hz)');
          
    %Filtro Cheby2
    Rs = 50; 
    H2 = FiltroCheby2(Fs, fc, NFFT, Rs, n, 'stopp', 0.15, 0.85);
    S_filt2 = H2.*S;
    s2 = ifft(S_filt2);
    s2 = s2(1:N)*(N/2);

    subplot(2,3,4);
    plot(t, s2);
    title('Sinal filtrado pelo Filtro Cheby2 rejeita banda');
    xlabel('Tempo (s)');
    ylabel('|S_filt| /mV');
    
##    S_f2 = fft(s2, NFFT)*2/N;
##
##    plot(f, abs(S_f2));
##    axis([0 Fs/2]);
##    title('Amplitude Espectral de S em função da Frequência');
##    ylabel('|S| /mV');
##    xlabel('Frequência (Hz)');
          

    %Filtro Ellip
    Rp1 = 5;
    Rs1 = 50;
    H3 = FiltroEllip(Fs, fc, NFFT, Rp1, Rs1, n, 'stopp', 0.15, 0.85);
    S_filt3 = H3.*S;
    s3 = ifft(S_filt3);
    s3 = s3(1:N)*(N/2);

    subplot(2,3,5);
    plot(t, s3);
    title('Sinal filtrado pelo Filtro ellip rejeita banda');
    xlabel('Tempo (s)');
    ylabel('|S_filt| /mV');
    
##    S_f3 = fft(s3, NFFT)*2/N;
##
##    plot(f, abs(S_f3));
##    axis([0 Fs/2]);
##    title('Amplitude Espectral de S em função da Frequência');
##    ylabel('|S| /mV');
##    xlabel('Frequência (Hz)');
          

    %Filtro FIR1
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
    ecg_fi = filtfilt(b,a,ecg);
    ##freqz(b,a);

    subplot(2,3,6);
    plot(t, ecg_fi);
    title('Filtro FIR Passa-banda 60 Hz');
    xlabel('Tempo /s');
    ylabel('Amplitude/mV');

##    S_fi = fft(ecg_fi, NFFT)*2/N;
## 
##    plot(f, abs(S_fi));
##    axis([0 Fs/2]);
##    title('Amplitude Espectral de S em função da Frequência');
##    ylabel('|S| /mV');
##    xlabel('Frequência (Hz)');

  case 2
      %Sinal original
      subplot(2,3,1);
      plot(t, val);
      title('Sinal, s(n), corrompido com ruído');
      xlabel('Tempo /s');
      ylabel('Amplitude /mV');
      
      %Filtro Butter
      H = FiltroButter(Fs, Fc, NFFT, n, 'altoo');
      S_filt = H.*S;
      s = ifft(S_filt);
      s = s(1:N)*(N/2);

      subplot(2,3,2);
      plot(t, s);
      title('Sinal filtrado pelo Filtro Butter passa alto');
      xlabel('Tempo (s)');
      ylabel('|S_filt| /mV');
      
##      S_f = fft(s, NFFT)*2/N;
##
##      plot(f, abs(S_f));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');

      %Filtro Cheby1
      Rp = 0.5;
      H1 = FiltroCheby1(Fs, Fc, NFFT, Rp, n, 'altoo');
      S_filt1 = H1.*S;

      s1 = ifft(S_filt1);
      s1 = s1(1:N)*(N/2);

      subplot(2,3,3);
      plot(t, s1);
      title('Sinal filtrado pelo Filtro Chevy1 passa alto');
      xlabel('Tempo (s)');
      ylabel('|S_filt| /mV');
      
##      S_f1 = fft(s1, NFFT)*2/N;
##
##      plot(f, abs(S_f1));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');
    
      %Filtro Cheby2
      Rs = 20; 
      H2 = FiltroCheby2(Fs, Fc, NFFT, Rs, n, 'altoo');
      S_filt2 = H2.*S;

      s2 = ifft(S_filt2);
      s2 = s2(1:N)*(N/2);

      subplot(2,3,4);
      plot(t, s2);
      title('Sinal filtrado pelo Filtro Chevy2 passa alto');
      xlabel('Tempo(s)');
      ylabel('|S_filt| /mV');
          
      S_f2 = fft(s2, NFFT)*2/N;

##      plot(f, abs(S_f2));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');

      %Filtro Ellip
      Rp1 = 5;
      Rs1 = 50;
      H3 = FiltroEllip(Fs, Fc, NFFT, Rp1, Rs1, n, 'altoo');
      S_filt3 = H3.*S;

      s3 = ifft(S_filt3);
      s3 = s3(1:N)*(N/2);

      subplot(2,3,5);
      plot(t, s3);
      title('Sinal filtrado pelo Filtro ellip passa alto');
      xlabel('Tempo (s)');
      ylabel('|S_filt| /mV');

      
      S_f3 = fft(s3, NFFT)*2/N;

##      plot(f, abs(S_f3));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');
          
      %Filtro FIR1
      ordem = 1000;
      a = 1;
      b = fir1(ordem,1/(Fs/2),'high');

      % filtrar sinal
      ecg_f = filtfilt(b,a,ecg);
      subplot(2,3,6);
      plot(t, ecg_f);
      title('Filtro FIR Passa-alto 1 Hz');
      xlabel('Tempo /s');
      ylabel('Amplitude /mV');
      
      S_fi = fft(ecg_f, NFFT)*2/N;


##      plot(f, abs(S_fi));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');

  case 3
      Fc1 = 65;
      %Sinal original
      subplot(2,3,1);
      plot(t, val);
      title('Sinal, s(n), corrompido com ruído');
      xlabel('Tempo /s');
      ylabel('Amplitude /mV');
      
      %Filtro Butter
      H = FiltroButter(Fs, Fc1, NFFT, n, 'baixo');
      S_filt = H.*S;
      s = ifft(S_filt);
      s = s(1:N)*(N/2);

      subplot(2,3,2);
      plot(t, s);
      title('Sinal filtrado pelo Filtro Butter passa baixo');
      xlabel('Tempo (s)');
      ylabel('|S_filt| /mV');

      
##      S_f = fft(s, NFFT)*2/N;
##
##      plot(f, abs(S_f));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');

      %Filtro Cheby1
      Rp = 0.5;
      H1 = FiltroCheby1(Fs, Fc1, NFFT, Rp, n, 'baixo');
      S_filt1 = H1.*S;

      s1 = ifft(S_filt1);
      s1 = s1(1:N)*(N/2);

      subplot(2,3,3);
      plot(t, s1);
      title('Sinal filtrado pelo Filtro Chevy1 passa baixo');
      xlabel('Tempo (s)');
      ylabel('|S_filt| /mV');

      S_f1 = fft(s1, NFFT)*2/N;

##      plot(f, abs(S_f1));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');

      %Filtro Cheby2
      Rs = 20; 
      H2 = FiltroCheby2(Fs, Fc1, NFFT, Rs, n, 'baixo');
      S_filt2 = H2.*S;

      s2 = ifft(S_filt2);
      s2 = s2(1:N)*(N/2);

      subplot(2,3,4);
      plot(t, s2);
      title('Sinal filtrado pelo Filtro Chevy2 passa baixo');
      xlabel('Tempo(s)');
      ylabel('|S_filt| /mV');
          
         
      S_f2 = fft(s2, NFFT)*2/N;

##      plot(f, abs(S_f2));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');
      
      
      %Filtro Ellip
      Rp1 = 5;
      Rs1 = 50;
      H3 = FiltroEllip(Fs, Fc1, NFFT, Rp1, Rs1, n, 'baixo');
      S_filt3 = H3.*S;

      s3 = ifft(S_filt3);
      s3 = s3(1:N)*(N/2);

      subplot(2,3,5);
      plot(t, s3);
      title('Sinal filtrado pelo Filtro ellip passa baixo');
      xlabel('Tempo (s)');
      ylabel('|S_filt| /mV');
      
      S_f3 = fft(s3, NFFT)*2/N;

##      plot(f, abs(S_f3));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');
          
      %Filtro FIR1
      ordem = 1000;
      a = 1;
      b = fir1(ordem,Fc1/(Fs/2),'low');

      % filtrar sinal
      ecg_f = filtfilt(b,a,ecg);
      subplot(2,3,6);
      plot(t, ecg_f);
      title('Filtro FIR Passa-baixo 65 Hz');
      xlabel('Tempo /s');
      ylabel('Amplitude /mV');
      
##      S_fi = fft(ecg_f, NFFT)*2/N;
##
##      plot(f, abs(S_fi));
##      axis([0 Fs/2]);
##      title('Amplitude Espectral de S em função da Frequência');
##      ylabel('|S| /mV');
##      xlabel('Frequência (Hz)');

 
endswitch
  

%--------------------------------------------------------------------------------------
