function H_Cheby1 = FiltroCheby1(fs, fc, NFFT, Rp, N, tipo, inf, sup)
  
  wp = fc/(fs/2);
  w = 0 : 2*pi/NFFT: 2*pi - 2*pi/NFFT;
  z = exp(1i*w);
  num = 1; 
  den = 1; 
  
  if (tipo == 'baixo')
    [zero,p,k] = cheby1(N,Rp, wp);
    for i = 1:N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    end
    H_Cheby1 = (num./den).*k; 
  endif
  
  if (tipo == 'stopp')
    [zero,p,k] = cheby1(N,Rp,[inf sup], 'stop');
    for i = 1:2.*N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    endfor
    H_Cheby1 = (num./den).*k; 
  endif
  
  if (tipo == 'altoo')
    [zero,p,k] = cheby1(N, Rp,(fc./(fs./2)), 'high');
    for i = 1:N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    endfor
    H_Cheby1 = (num./den).*k; 
  endif
  
  if (tipo == 'passa')
    [zero,p,k] = cheby1(N, Rp,([inf sup]./(fs./2)));
    for i = 1:2.*N
      num = (z-zero(i)).*num;
      den = (z-p(i)).*den;
    endfor
    H_Cheby1 = (num./den).*k; 
  endif
end 