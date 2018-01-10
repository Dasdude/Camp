function rho = wdensity(S,T,P)
%WDENSITY  Returns the water density 
%
%  CALL:  rho = wdensity(S,T,P);
%
%   rho =  water density [kg/m^3]
%   S   = salinity       [psu]       (Default 35)
%   T   = temperature    [degrees C] (default 4)
%   P   = pressure       [db]        (default 0) 
%
% Example:
%  S = linspace(20,40)'; T = linspace(4,20)';
%  [S1 T1]=meshgrid(S,T);
%  sc = contour(S,T,wdensity(S1,T1)); clabel(sc)
%  xlabel('Salinity'),ylabel('Temperature')

% REFERENCES:
%  Fofonoff, P. and Millard, R.C. Jr (1983)
%  Algorithms for computation of fundamental properties of 
%  seawater, Unesco Tech. Pap. in Mar. Sci., No. 44, 53, pp 17-18
%
%  Millero, F.J & Poisson, A. (1981)
%  International one-atmosphere equation of state for seawater.
%  Deep-Sea Research, Vol. 28A, No.6, pp 625-629. 

%"An overview of the 1995 SWARM shallow water internal
%	acoustic scattering experiment", Apel, et el. IEEE J.
%	Ocean. Eng., Vol 22, No. 3, July 1997. 

%	"A Surface-Trapped Intrusion of Slope Water onto the Continental 
%	Shelf in the Middle Atlantic Bight", Gawarkiewicz, et. el. ,
%	Geophys. Res. Let. 23(25), 3763-3766, 1996.

%rho = 1026; Old value

error(nargchk(0,3,nargin))
if nargin<1|isempty(S), S = 35; end
if nargin<2|isempty(T), T = 4; end
if nargin<3|isempty(P), P = 0; end
[ec, S,T,P] = comnsize(S,T,P);
if ec>0,
  error('S, T and P must be scalar or of common size')
end

% Density of fresh water at atmospheric pressure
%-------------------------------------------------
rho0 = 999.842594; 
rho = rho0(ones(size(S)));

k = find(T~=0);
if any(k),
  Tk = T(k);
  a  = [6.536332e-9, -1.120083e-6, 1.001685e-4, -9.095290e-3, 6.793952e-2];
  rho(k) = rho(k)+(a(5)+(a(4)+(a(3)+(a(2)+a(1)*Tk).*Tk).*Tk).*Tk).*Tk;
  %rho(k) = rho(k) + polyval(a,Tk);
end

% Density of salt water at atmospheric pressure:
%------------------------------------------------
k1 = find(S~=0);
if any(k1),
  c3 = -5.72466e-3;
  b5 = 8.24493e-1;
  tmp  = b5(ones(size(k1)));
  tmp1 = c3(ones(size(k1)));
  Tk   = T(k1);
  k12  = find(Tk~=0);
  if any(k12),
    Tk = Tk(k12);
    b = [ 5.3875e-9, -8.2467e-7,  7.6438e-5, -4.0899e-3, 0];
    c = [ -1.6546e-6, 1.0227e-4, 0];
    tmp(k12)  = tmp(k12)+(b(4) + (b(3) + (b(2) + b(1)*Tk).*Tk).*Tk).*Tk; 
    tmp1(k12) = tmp1(k12)+(c(2) + c(1)*Tk).*Tk;
    clear Tk
  end
  d = 4.8314e-4;
  Sk = S(k1);
  rho(k1) = rho(k1) + (tmp + tmp1.*sqrt(Sk) + d*Sk).*Sk;	       
end

k2 = find(P~=0);
if any(k2),
  Pk  = P(k2);
  K   = seck(S(k2),T(k2),Pk);
  Pk  = Pk/10;  % convert from db to atm pressure units
  rho(k2) = rho(k2)./(1-Pk./K);
end
return


function K = seck(S,T,P)
%SECK  Secant bulk modulus (K) of sea water
%
% CALL:  K = seck(S,T,P);
%
%   S = salinity    [psu]
%   T = temperature [degrees C]
%   P = pressure    [db]
%
%   K = Secant Bulk Modulus  [bars]
% 
%    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980. 
%    UNESCO polynomial implementation.


error(nargchk(3,3,nargin))


% Fresh water terms of K at atmosphere pressure.
%------------------------------------------------

K0 = 19652.21;
K  = K0(ones(size(S)));

k1 = find(T~=0);
if any(k1),
  Tk = T(k1);
  e = [-5.155288E-5, 1.360477E-2, -2.327105, 148.4206,0];
  K(k1)  = K(k1) + (e(4) + (e(3) + (e(2) + e(1)*Tk).*Tk).*Tk).*Tk;   % eqn 19
end



% Salt water terms of K at atmosphere pressure
%---------------------------------------------

SR = sqrt(S);
k1 = find(S~=0);
if any(k1),
  Tk = T(k1);
  f = [-6.1670E-5, 1.09987E-2,-0.603459,54.6746];
  g = [-5.3009E-4, 1.6483E-2,7.944E-2];
  
  K(k1) = K(k1) + (  f(4) + (f(3) + (f(2) + f(1)*Tk).*Tk).*Tk ...
      + (g(3) + (g(2) + g(1)*Tk).*Tk).*SR(k1)).*S(k1);      % eqn 16
end


k2 = find(P~=0);
if any(k2),
  A0 = 3.239908; B0 = 8.50935E-5;
  A  = A0(ones(size(k2)));
  B  = B0(ones(size(k2)));
  Tk = T(k2);
  
  k12 = find(Tk~=0);
  if any(k12)
    Tk1 = Tk(k12);
    k12 = k(k12);
    h  = [-5.77905E-7, 1.16092E-4, 1.43713E-3, 0];
    A(k12)  = A(k12) + (h(3) + (h(2) + h(1).*Tk1).*Tk1).*Tk1;
    
    b = [ 5.2787E-8, -6.12293E-6,0];
    B(k12)  = B(k12) + (b(2) + b(1)*Tk1).*Tk1;
  end
  J = 1.91075E-4;
  I =[ -1.6078E-6, -1.0981E-5  2.2838E-3];
  A  = A + (I(3) + (I(2) + I(1)*Tk).*Tk + J*SR(k2)).*S(k2); 
  
  m = [ 9.1697E-10, 2.0816E-8, -9.9348E-7];
  B = B + (m(3) + (m(2) + m(1)*Tk).*Tk).*S(k2);  
  
  Pk = P(k2)/10;  %convert from db to atmospheric pressure units
  K(k2) = K(k2) + (A + B.*Pk).*Pk;  % eqn 15
end
return









