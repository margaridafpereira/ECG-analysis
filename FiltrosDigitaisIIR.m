%amplitude espectral 1

clear all
close all
clc
%-------------------------------------------------------
pkg load signal
%%------------------------------------------------------
%% inicializa constantes e variaveis

f1 = 12.5; %frequencia de s1, em Hz (componente fundamental)
A1 = 5; %amplitude de s1, em V
fi1 = 0; %fase inicial de s1

f2 = 150; %frequencia de s2, em Hz
A2 = 5; %amplitude de s2, em V
fi2 = pi/2; %fase inicial de s2

f3 = 750; %frequencia de s3, em Hz
A3 = 5; %amplitude de s3, em V
fi3 = pi/4; %fase inicial de s3

T1 = 1/f1; %periodo da componente s1
Nciclos = 2.5; %numero de ciclos completos que s1 deve descrever

t_total = Nciclos * T1; 
fs = 2000; %frequencia de amostragem
n = t_total*fs; %numero de amostras

fc = 100; %frequencia de corte, em Hz
%%---------------------------------------------------------
%% calcula componentes do sinal
ts = t_total/n; %periodo de amostragem 
t = 0: ts: (t_total-ts); %tempo

w1 = 2*pi*f1; %%frequencia angular para s1
s1 = A1 * sin(w1*t + fi1); 

w2 = 2*pi*f2; %frequencia angular para s2
s2 = A2 * sin(w2*t + fi2);

w3 = 2*pi*f3; %frequencia angular para s3
s3 = A3 * sin(w3*t + fi3); 

%ruido branco, criado com a função rand
Ar = 1; 
rbr = 2.*Ar.*rand(1,n) - Ar;
##plot(t, rbr);
##title('Ruido Branco');
##xlabel('Tempo /s');
##ylabel('Amplitude /mV');
s = s1 + s2 + s3 + rbr; #s é a soma de s1, s2 e s3
%-----------------------------------------------------------------
################################
#####REPRESENTAÇÃO DO SINAL#####
################################

##plot(t, s); #reprsentação do sinal em função do tempo
##title('Sinal, s(n), corrompido com ruído');
##xlabel('Tempo /s');
##ylabel('Amplitude /mV');

NFFT=2^nextpow2(n); %próxima potência de 2 do comprimento de s
f = 0: fs/NFFT: fs- fs/NFFT; 
S = fft(s, NFFT)*2/n; 

w = 0 : 2*pi/NFFT: 2*pi - 2*pi/NFFT;

%------------------------------------------------------------------------------------
##############################
#####Filtros digitais IIR#####
##############################
%------------------------------------------------------------------------
N = 4;#ordem do filtro 

choice = menu('Sinal', 'Original', 'Filtrado')
switch choice
  case 1
    
    %representação no dominio das frequencias
    plot(f, abs(S)); 
    axis ([0 fs/2]);
    title('Amplitude Espectral de S em função da Frequência');
    xlabel('Frequência (Hz)');
    ylabel('|S| /mV');
    
  case 2
  choice1 = menu('Tipo de Filtro', 'Butter', 'Cheby1', 'Cheby2', 'ellip')
  switch choice1
    case 1
      choice2 = menu('Tipo de Filtro', 'Passa-baixo', 'Stop', 'Passa-Alto', 'Passa-Banda')
      switch choice2
        case 1
          H = FiltroButter(fs, fc, NFFT, N, 'baixo');
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Butter passa baixo');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
          
       case 2
          H = FiltroButter(fs, fc, NFFT, N, 'stopp', 0.15, 0.85);
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Butter rejeita banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
      case 3
          H = FiltroButter(fs, fc, NFFT, N, 'altoo');
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Butter passa alto');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
      case 4
          H = FiltroButter(fs, fc, NFFT, N, 'passa', 500, 1000);
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Butter passa banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
      endswitch
    case 2
      choice3 = menu('Tipo de Filtro', 'Passa-baixo', 'Stop', 'Passa-Alto', 'Passa-Banda')
      switch choice3
        case 1
          Rp = 10;
          H = FiltroCheby1(fs, fc, NFFT, Rp, N, 'baixo');
          S_filt = H.*S;
          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro cheby1 passa baixo');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
      case 2
          Rp = 5; 
          H = FiltroCheby1(fs, fc, NFFT, Rp, N, 'stopp', 0.15, 0.85);
          S_filt = H.*S;
          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Cheby1 rejeita banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
      case 3
          Rp = 0.5;
          H = FiltroCheby1(fs, fc, NFFT, Rp, N, 'altoo');
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Chevy1 passa alto');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
          
      case 4
          Rp = 3;
          H = FiltroCheby1(fs, fc, NFFT,Rp, N, 'passa', 500, 1000);
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Cheby1 passa banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
        
        endswitch
      case 3
        choice4 = menu('Tipo de Filtro', 'Passa-baixo', 'Stop', 'Passa-Alto', 'Passa-Banda')
        switch choice4
          case 1
          Rs = 40; 
          H = FiltroCheby2(fs, fc, NFFT, Rs, N, 'baixo');
          S_filt = H.*S;
          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro cheby2 passa baixo');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
        case 2
          Rs = 50; 
          H = FiltroCheby2(fs, fc, NFFT, Rs, N, 'stopp', 0.15, 0.85);
          S_filt = H.*S;
          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Cheby2 rejeita banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
        case 3
          Rs = 20; 
          H = FiltroCheby2(fs, fc, NFFT, Rs, N, 'altoo');
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Chevy2 passa alto');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
          
        case 4
          Rs = 40;
          H = FiltroCheby2(fs, fc, NFFT,Rs, N, 'passa', 500, 1000);
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Cheby2 passa banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
          
        endswitch
      case 4
        choice5 = menu('Tipo de Filtro', 'Passa-baixo', 'Stop', 'Passa-Alto', 'Passa-Banda')
        switch choice5
          case 1
          Rp = 5;
          Rs = 40;
          H = FiltroEllip(fs, fc, NFFT, Rp, Rs, N, 'baixo');
          S_filt = H.*S;
          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro ellip passa baixo');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
        case 2
          Rp = 5;
          Rs = 50;
          H = FiltroEllip(fs, fc, NFFT, Rp, Rs, N, 'stopp', 0.15, 0.85);
          S_filt = H.*S;
          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro ellip rejeita banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
        case 3
          Rp = 3;
          Rs = 50;
          H = FiltroEllip(fs, fc, NFFT, Rp, Rs, N, 'altoo');
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro ellip passa alto');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
          
        case 4
          Rp = 3;
          Rs = 40;
          H = FiltroEllip(fs, fc, NFFT, Rp, Rs, N, 'passa', 500, 1000);
          S_filt = H.*S;

          plot(f, abs(S), f, abs(S_filt));
          title('Sinal Original e sinal filtrado pelo Filtro Ellip passa banda');
          xlabel('Frequência (Hz)');
          ylabel('|S_filt| /mV');
          legend('Sinal Original', 'Sinal Filtrado');
         endswitch
      endswitch
    
endswitch
  

%--------------------------------------------------------------------------------------
