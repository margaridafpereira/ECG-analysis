function H_Ellip = FiltroEllip(fs, fc, NFFT, Rp1, Rs1, N, tipo, inf, sup)
  
  wp = fc/(fs/2);
  w = 0 : 2*pi/NFFT: 2*pi - 2*pi/NFFT;
  z = exp(1i*w);
  num = 1; 
  den = 1; 
  
  if (tipo == 'baixo')
    [zero,p,k] = ellip(N,Rp1,Rs1, wp);
    for i = 1:N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    end
    H_Ellip = (num./den).*k; 
  endif
  
  if (tipo == 'stopp')
    [zero,p,k] = ellip(N,Rp1,Rs1,[inf sup], 'stop');
    for i = 1:2.*N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    endfor
    H_Ellip = (num./den).*k; 
  endif
  
  if (tipo == 'altoo')
    [zero,p,k] = ellip(N, Rp1, Rs1,(fc./(fs./2)), 'high');
    for i = 1:N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    endfor
    H_Ellip = (num./den).*k; 
  endif
  
  if (tipo == 'passa')
    [zero,p,k] = ellip(N, Rp1, Rs1,([inf sup]./(fs./2)));
    for i = 1:2.*N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    endfor
    H_Ellip = (num./den).*k; 
  endif
end 