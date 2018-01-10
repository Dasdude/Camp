 function [int,tol] = ccquad(fun,a,b,n,trace)
% CCQUAD Numerical integration using a Clenshaw-Curtis quadrature.
%
% CALL:	 [int, tol] = ccquad(Fun,a,b,n) 
%
%	int = evaluated integral
%       tol = estimate of the absolute error (usually conservative).
%	Fun = string containing the name of the function or a directly 
%             given expression enclosed in parenthesis. 
%       a,b = integration limits
%         n = number of base points (abscissas). Default n=10
%
%  The integral is exact for polynomials of degree n or less.
%  Usually this routine gives accurate answers for smooth functions.  
%
%Example:% Integration of exp(x) from 0 to 1:
%
%   a=0; b=1;
%   ccquad('exp(x)',a,b) % Should be exp(1)-1
%
% See also  gaussq, qrule2d

% References:
% [1] Goodwin, E.T. (1961),
% "Modern Computing Methods",
% 2nd edition, New yourk: Philosophical Library, pp. 78--79
%
% [2] Clenshaw, C.W. and Curtis, A.R. (1960),
% Numerische Matematik, Vol. 2, pp. 197--205

% tested on: matlab 5.3
% history:
% revised pab 22Nov2004
% Added the possibility of using a function handle.
% by	Per A. Brodtkorb 26.07.1999

if nargin<5|isempty(trace),
   trace=0;
end
if nargin<4|isempty(n),
   n = 10;
else
   % make sure n is even
   n = 2*ceil(n/2);
end;

if (isa(fun,'char') &  any(fun=='(')), %  & any(fun=='x'),
  exec_string=['f=',fun ';']; %the call function is already setup
else
  if isa(fun,'function_handle')
    fun = func2str(fun);
  end
  %setup string to call the function
  exec_string=['f=feval(fun,x);'];
end




s = (0:n)';
s2 =(0:2:n)';

x = cos(pi*s/n)*(b-a)/2+(b+a)/2;
eval(exec_string);
if trace==1,
  plot(x,f,'+')
end

% using a Gauss-Lobatto variant, i.e., first and last
% term f(a) and f(b) is multiplied with 0.5
f(1) = f(1)/2;
f(n+1) = f(n+1)/2;

if 1,%fft for faster calculations
 c=real(fft(f(1:n)));
 c=2/n*(c(1:n/2+1)+f(n+1)*cos(pi*s2));
else %old call: slow for large n
  c = 2/n * cos(s2*s'*pi/n) * f;
end
c(1) = c(1)/2;
c(n/2+1) = c(n/2+1)/2;

int = (a-b)*sum(c./(s2-1)./(s2+1) );
tol = abs(c(n/2+1));

