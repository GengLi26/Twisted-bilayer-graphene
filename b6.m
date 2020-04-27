function result = b6(kx,ky,l,m,a_moire)
  
%% calculate phi_lm 
  % l, m  numbers; 
  % a_moire is moire lattice constant

  result=exp(-i*a_moire*(l*(sqrt(3)/2)*kx+(m-l/2)*ky));

end