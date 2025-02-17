function h=hldpi2fft(A,kernel,L)
%HLDPI2FFT L-stage DPI estimate of smoothing parameter for 2D data
%           (DPI=Direct Plug-In)
%
% CALL  h=hldpi2fft(A,kernel,L) 
%  
%       hs = smoothing parameter
%   data   = data matrix, size N x D (D = # dimensions ==2)
%   kernel = 'gaussian'      - Gaussian kernel
%            'epanechnikov'  - Epanechnikov kernel.
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangluar'    - Triangular kernel.
%            'rectangular'   - Rectanguler kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.  
%        L = 0,1,2,3   (default 2)
%  
%  Note that only the first 4 letters of the kernel name is needed.
%
%  Example: 
%   x  = wnormrnd(0,1,50,2);
%   hs = hldpi2fft(x,'gauss',1);
%
% See also  hste, hbcv, hboot, hos, hlscv, hscv, hstt, kde, kdefun

%  Wand, M.P. and Jones, M.C. (1994) 
%  "Multivariate Plug-In Bandwidth Selection", 
%  Computational statistics

% tested on: matlab 5.2
% history:
% revised pab aug2005
% -added ad hoc support for other kernels than Gaussian
% revised pab dec2003
%  -fixed some bugs  
%  -added todo comments 
%  -added call to gridcount  
% -updated example  
% -fixed a bug  
% revised pab 3.11.1999
% from kdetools  by  Christian C. Beardah 1996




[n,d]=size(A);

if d~=2,
  error('Wrong size of data')
end
if nargin<2 |isempty(kernel)
  kernel = 'gauss';
end

if nargin<3|isempty(L),
  L=2;
else
  L=min(abs(L),3);
end;

if any(strcmpi(kernel(1:4),{'gaus','norm'})),
  r = 4;
else,
    r=1;
end;

inc=128;

P=2*inc;

xmin = min(A,[],1);      % Find the minimum value of x.
xmax = max(A,[],1);      % Find the maximum value of x.
xrange=xmax-xmin; % Find the range of x.

% xa holds the x 'axis' vector, defining a grid of x values where 
% the k.d. function will be evaluated and plotted.

ax=xmin-xrange/4;
bx=xmax+xrange/4;

deltax = (bx-ax)/(inc-1);

X =[linspace(ax(1),bx(1),inc).' , linspace(ax(2),bx(2),inc).'];

% Find the binned kernel weights, c.
c = gridcount(A,X);

% R= int(mkernel(x)^2)
% mu2= int(x^2*mkernel(x))
kernel2 = 'gauss';
% values for gaussian kernel
[mu2ns,Rns] = kernelstats(kernel2);
Rns = Rns^d;
STEconstant2 = Rns /(mu2ns^(2)*n);

[mu2,R] = kernelstats(kernel);
R = R^d;
STEconstant = R /(mu2^(2)*n);




  K22=1/(2*pi);
  K40=3/(2*pi);
  K04=3/(2*pi);
  
  K42=-3/(2*pi);
  K24=-3/(2*pi);
  K60=-15/(2*pi);
  K06=-15/(2*pi);

  K44=9/(2*pi);
  K62=15/(2*pi);
  K26=15/(2*pi);
  K80=105/(2*pi);
  K08=105/(2*pi);

sigma1=std(A(:,1));
sigma2=std(A(:,2));

rho = corrcoef(A);
rho = rho(1,2);

if L==3,

  lam100=945;
  lam010=945;
  lam82=105+840*rho^2;
  lam28=105+840*rho^2;
  lam64=45+540*rho^2+360*rho^4;
  lam46=45+540*rho^2+360*rho^4;

  psi100=lam100/(128*pi*sigma1^11*sigma2*(1-rho^2)^5.5);
  psi010=lam010/(128*pi*sigma1*sigma2^11*(1-rho^2)^5.5);
  psi82=lam82/(128*pi*sigma1^9*sigma2^3*(1-rho^2)^5.5);
  psi28=lam28/(128*pi*sigma1^3*sigma2^9*(1-rho^2)^5.5);
  psi64=lam64/(128*pi*sigma1^7*sigma2^5*(1-rho^2)^5.5);
  psi46=lam46/(128*pi*sigma1^5*sigma2^7*(1-rho^2)^5.5);

% The following are the smoothing parameters for 
% estimation of psi80, psi08, psi62, psi26 and psi44:

  a80=(-2*K80/(mu2ns*(psi100+psi82)*n))^(1/12);
  a08=(-2*K08/(mu2ns*(psi28+psi010)*n))^(1/12);
  a62=(-2*K62/(mu2ns*(psi82+psi64)*n))^(1/12);
  a26=(-2*K26/(mu2ns*(psi46+psi28)*n))^(1/12);
  a44=(-2*K44/(mu2ns*(psi64+psi46)*n))^(1/12);

  psi80=psim('80',A,a80);
  psi08=psim('08',A,a08);
  psi62=psim('62',A,a62);
  psi26=psim('26',A,a26);
  psi44=psim('44',A,a44);

end;

if L==2 | L==3,

  if L==2,
    lam80=105;
    lam08=105;
    lam62=15+90*rho^2;
    lam26=15+90*rho^2;
    lam44=9+72*rho^2+24*rho^4;

    psi80=lam80/(64*pi*sigma1^9*sigma2*(1-rho^2)^4.5);
    psi08=lam08/(64*pi*sigma1*sigma2^9*(1-rho^2)^4.5);
    psi62=lam62/(64*pi*sigma1^7*sigma2^3*(1-rho^2)^4.5);
    psi26=lam26/(64*pi*sigma1^3*sigma2^7*(1-rho^2)^4.5);
    psi44=lam44/(64*pi*sigma1^5*sigma2^5*(1-rho^2)^4.5);
  end;

% The following are the smoothing parameters for 
% estimation of psi60, psi06, psi42 and psi24:

  a60=(-2*K60/(mu2ns*(psi80+psi62)*n))^0.1;
  a06=(-2*K06/(mu2ns*(psi26+psi08)*n))^0.1;
  a42=(-2*K42/(mu2ns*(psi62+psi44)*n))^0.1;
  a24=(-2*K24/(mu2ns*(psi44+psi26)*n))^0.1;

  psi60=psim('60',A,a60);
  psi06=psim('06',A,a06);
  psi42=psim('42',A,a42);
  psi24=psim('24',A,a24);

end;

if L==1 | L==2 | L==3,

  if L==1,
    lam06=15;
    lam60=15;
    lam24=3+12*rho^2;
    lam42=3+12*rho^2;

    psi60=-lam60/(32*pi*sigma1^7*sigma2*(1-rho^2)^3.5);
    psi06=-lam06/(32*pi*sigma1*sigma2^7*(1-rho^2)^3.5);
    psi42=-lam42/(32*pi*sigma1^5*sigma2^3*(1-rho^2)^3.5);
    psi24=-lam24/(32*pi*sigma1^3*sigma2^5*(1-rho^2)^3.5);
  end;

% The following are the smoothing parameters for 
% estimation of psi40, psi04 and psi22:

  a40=(-2*K40/(mu2ns*(psi60+psi42)*n))^0.125;
  a04=(-2*K04/(mu2ns*(psi24+psi06)*n))^0.125;
  a22=(-2*K22/(mu2ns*(psi42+psi24)*n))^0.125;

% Brute force technique:

  psi40=psim('40',A,a40);
  psi04=psim('04',A,a04);
  psi22=psim('22',A,a22);

% Fourier technique:

% Psi40:

  h=real(a40);

% Obtain the kernel weights.

  L1=max(floor(r*h./deltax));

  L2=min(L1,inc-1);

  kw=zeros(2*inc,2*inc);

  [X1,Y1]=meshgrid(deltax(1)*[-L2:L2]/h,deltax(2)*[-L2:L2]/h);
  idx(1:d) = {(inc-L2+1):(inc+L2+1)};
  
  kw(idx{:}) = deriv2(X1,Y1,'40')/(n*h^2);

% Apply 'fftshift' to kw.

  kw = fftshift(kw);

% Perform the convolution.

  z = real(ifft2(fft2(c,P,P).*fft2(kw)));

  z=z(1:inc,1:inc).*(z(1:inc,1:inc)>0);

  psi40=sum(sum(c.*z(1:inc,1:inc)))/(n^2*h^5);

end;

if L==0,

  lam04=3;
  lam40=3;
  lam22=3/2;

  psi40=-lam40/(16*pi*sigma1^5*sigma2*(1-rho^2)^2.5);
  psi04=-lam04/(16*pi*sigma1*sigma2^5*(1-rho^2)^2.5);
  psi22=-lam22/(16*pi*sigma1^3*sigma2^3*(1-rho^2)^2.5);

end;



h1=(psi04^0.75*R/(n*mu2^2*psi40^0.75*(sqrt(psi40*psi04)+psi22)))^(1/6);
h2=h1*(psi40/psi04)^0.25;

h=real([h1, h2]);


return


function p=psim(m,A,h)

% PSIM        Calculate psi_m by brute force.
%
%             PSIM(M,A,H)
%
%             Christian C. Beardah 1996

n=length(A);

Xi=zeros(n,n);
Xj=Xi;
Yi=Xi;
Yj=Xi;

o=ones(1,n);

Xj=A(:,1)*o;
Xi=Xj';
Yj=A(:,2)*o;
Yi=Yj';

p=sum(sum(deriv2((Xi-Xj)/h,(Yi-Yj)/h,m) ))/(h*n^2);
return
