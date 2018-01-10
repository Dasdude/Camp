function [int, tol1,ix]= gaussq(fun,xlow,xhigh,tol,trace,varargin)
%GAUSSQ Numerically evaluates a integral using a Gauss quadrature.
%       The Quadrature integrates a (2m-1)th order  polynomial exactly 
%       and the  integral is of the form
%              b                              
%             Int (p(x)* Fun(x)) dx 
%              a                 
%
% CALL:
% [int, tol] = gaussq('Fun',xlow,xhigh,[reltol wfun],[trace,gn],p1,p2,....)
% [int, tol] = gaussq('Fun',xlow,xhigh,[reltol wfun],[trace,gn],alpha,p1,p2,....)
% [int, tol] = gaussq('Fun',xlow,xhigh,[reltol wfun],[trace,gn],alpha,beta,p1,p2,....)
%
%       int = evaluated integral
%       tol = absolute tolerance abs(int-intold)
%       Fun = string containing the name of the function or a directly given 
%             expression enclosed in parenthesis. The function may depend on
%             parameters alpha,beta, pi.
%xlow,xhigh = integration limits
%    reltol = relative tolerance default=1e-3
%      wfun = weight function
%         1  p(x)=1                       a =-1,   b = 1 Legendre (default)
%         2  p(x)=1/sqrt((x-a)*(b-x)),    a =-1,   b = 1 Chebyshev of the
%                                                             first kind
%         3  p(x)=sqrt((x-a)*(b-x)),      a =-1,   b = 1 Chebyshev of the 
%                                                            second kind
%         4  p(x)=sqrt((x-a)/(b-x)),      a = 0,   b = 1
%         5  p(x)=1/sqrt(b-x),            a = 0,   b = 1
%         6  p(x)=sqrt(b-x),              a = 0,   b = 1
%         7  p(x)=(x-a)^alpha*(b-x)^beta  a =-1,   b = 1 Jacobi 
%                                     alpha, beta >-1 (default alpha=beta=0)
%         8  p(x)=x^alpha*exp(-x)         a = 0,   b = inf generalized Laguerre
%         9  p(x)=exp(-x^2)               a =-inf, b = inf Hermite
%        10  p(x)=1                       a =-1,   b = 1 Legendre (slower than 1)
%
%   trace = for non-zero TRACE traces the function evaluations 
%              with a point plot of the integrand.
%      gn = number of base points and weight points to start the 
%            integration with (default=2)
%p1,p2,...= coefficients to be passed directly to function Fun:   
%                  G = Fun(x,p1,p2,...).
%
% Note that int is the common size of xlow, xhigh and p1,p2,....
% Example:% a) integration of x.^2 from 0 to 2 and from 1 to 4
%         % b) integration of x^2*exp(-x) from zero to infinity 
%
%    gaussq('(x.^2)',[0 1],[2 4])             % a)
%    gaussq('(1)',0,inf,[1e-3 8],[],2)        % b)
%    gaussq('(x.^2)',0,inf,[1e-3 8],[],0)     % b)
%
% See also  qrule, gaussq2d


% References 
%
% [1]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [2] Davis and Rabinowitz (1975) 'Methods of Numerical Integration', page 365,
%     Academic Press.
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% [4] Abromowitz and Stegun (1954) 'Handbook of mathematical functions'




% tested on: Matlab 5.3
% history:
% Revised pab 22Nov2004
% -Added the possibility of using a function handle.  
% Revised pab 09.09.2002
% -added the possibility of using a inline function.
% revised pab 27.03.2000
%  - fixed a bug for p1,p2,... changed to varargin in input
% revised pab 19.09.1999
%   documentation
%  by Per A. Brodtkorb 30.03.99, 17.02.99 :
%   wfun 1: from grule.m in NIT toolbox, see ref [2] 
%   wfun 2-6: see ref [4]
%   wfun 7-10:  Adapted from Netlib routine gaussq.f see ref [1,3]
%  -accept multiple integrationlimits, int is the common size of xlow and xhigh
%  -optimized by only computing the integrals which did not converge.
%  -enabled the integration of directly given functions enclosed in 
%     parenthesis. Example: integration from 0 to 2 and from 2 to 4 for x is done by:
%                        gaussq('(x.^2)',[0 2],[2 4])

global ALPHA1 BETA1 ALPHA2
wfun=1;
if nargin<4| isempty(tol),
  tol=1e-3;
elseif length(tol)==2,
  wfun=tol(2);
   tol=tol(1);
end

P0=varargin;
NP=length(P0);

istart      = 0; 
alpha1      = 0;
beta1       = 0;
FindWeights = 1;

x = [];
y = [];

switch wfun
  case 7,
    istart=2;
    if ((NP>=1) & (~isempty(P0{1}))), alpha1 = P0{1};    end
    if ((NP>=2) & (~isempty(P0{1}))), beta1  = P0{2};    end
    if isempty(ALPHA1)|isempty(BETA1),
    elseif ALPHA1==alpha1 & BETA1==beta1,
      FindWeights=0;
    end
    ALPHA1=alpha1;BETA1=beta1;
    %remember what type of weights are saved as global constants
  case 8,
    istart=1;
    if ((NP>=1) & (~isempty(P0{1}))), alpha1 = P0{1};    end
    if isempty(ALPHA2),
    elseif ALPHA2==alpha1,
      FindWeights=0;
    end
    ALPHA2=alpha1;
    %remember what type of weights are saved as global constants
  otherwise,FindWeights=0;
end

P0(1:istart)=[];

gn=2;
if nargin <5|isempty(trace) ,
  trace = 0; 
elseif length(trace)==2,
  gn    = trace(2);
  trace = trace(1);
end


if prod(size(xlow))==1,% make sure the integration limits have correct size
  xlow=xlow(ones(size(xhigh)));;
elseif prod(size(xhigh))==1,
  xhigh=xhigh(ones(size(xlow)));;
elseif any( size(xhigh)~=size(xlow) )
  error('The input must have equal size!')
end
[N M]=size(xlow);%remember the size of input

num_parameters=NP-istart;
if num_parameters>0,
  ptxt = sprintf('p%d,',1:num_parameters);
  ptxt(end)=[]; % remove ','
  eval(sprintf('[%s]=deal(P0{:});',ptxt));
end
%Old call
%for ix=1:num_parameters,
%  eval(['p' int2str(ix) '=P0{ix};']); % variable # i
%end

if (isa(fun,'char') &  any(fun=='(')), %  & any(fun=='x'),
  exec_string=['y=',fun ';']; %the call function is already setup
else
  %setup string to call the function
  if isa(fun,'function_handle')
    fun = func2str(fun);
  end
  exec_string=['y=feval(fun,x'];
  for ix=1:num_parameters,
    xvar=['p' int2str(ix)]; % variable # i
    if eval(['isnumeric(' ,xvar,')  & length(',xvar,'(:)) ~=1' ]) ,
	if N*M==1,     
	  eval(['[N M]=size(', xvar, ');']); 
	elseif  eval(['N*M~=prod(size(', xvar, '))']),
	  error('The input must have equal size!')
	end	
      eval([xvar, '=' ,xvar ,'(:);']); %make sure it is a column 
      exec_string=[exec_string,',' xvar '(k,ones(1,gn) )']; %enable integration with multiple arguments
    else
      exec_string=[exec_string,',' xvar];
    end
  end
  exec_string=[exec_string,');'];
end


nk   = N*M;% # of integrals we have to compute
k    = (1:nk)';
int  = zeros(nk,1);
tol1 = int;

wtxt    = sprintf('%d_%d',gn,wfun); % number of weights and
                                    % weightfunction used
cbtxt = sprintf('cb%s',wtxt); %base points string
cwtxt = sprintf('cw%s',wtxt); %weights string
eval(sprintf('global %s %s ;',cbtxt,cwtxt));
if isempty(eval(['cb',wtxt]))|FindWeights ,  
  % calculate the weights if necessary
  eval(sprintf('[%s,%s]=qrule(gn,wfun,alpha1,beta1);',cbtxt,cwtxt));
end

%setup mapping parameters and execution strings
xlow=xlow(:);
jacob=(xhigh(:)-xlow(:))/2;

switch wfun,% this is clumsy and can written more tidy
  case {1 ,10},
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k);'];
  case 2,  
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2);'];
  case 3, 
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k).^2;'];
  case 4,  
    calcx_string=['(ones(nk,1),:)).*jacob(k,ones(1, gn ))*2+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k)*2;'];
  case 5,  
    calcx_string=['(ones(nk,1),:)).*jacob(k,ones(1, gn ))*2+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*sqrt(jacob(k)*2);'];
  case 6,  
    calcx_string=['(ones(nk,1),:)).*jacob(k,ones(1, gn ))*2+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*sqrt(jacob(k)*2).^3;'];
  case 7, 
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k).^(alpha1+beta1+1);'];
  case {8,9}
     calcx_string=['(ones(nk,1),:));'];
     int_string=['(ones(nk,1),:).*y,2);'];
otherwise error('unknown option')
end


eval(['x=(',cbtxt, calcx_string]); % calculate the x values            
eval(exec_string); % calculate function values  y=fun(x)                        
eval(['int(k)=sum(',cwtxt, int_string]);       % do the integration
int_old=int;

if trace==1,
  x_trace=x(:);
  y_trace=y(:);
end

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to int  and  compared to tol)
%  2) There are more than 11 iterations. This should NEVER happen. 

converge='n';
for i=1:10,
  gn=gn*2;% double the # of weights
  wtxt = sprintf('%d_%d',gn,wfun); % # of weights and weight function used
  %eval(['global cb',wtxt,' cw',wtxt]);
  cbtxt = sprintf('cb%s',wtxt); %base points string
  cwtxt = sprintf('cw%s',wtxt); %weights string
  eval(sprintf('global %s %s ;',cbtxt,cwtxt));
  if isempty(eval(cbtxt))|FindWeights ,  
    % calculate the weights if necessary
     eval(sprintf('[%s,%s]=qrule(gn,wfun,alpha1,beta1);',cbtxt,cwtxt));
  end
 
  eval(['x=(',cbtxt, calcx_string]); % calculate the x values  
  eval(exec_string);                 % calculate function values  y=fun(x)
  eval(['int(k)=sum(',cwtxt, int_string]);% do the integration
  
  if trace==1,
    x_trace=[x_trace;x(:)];
    y_trace=[y_trace;y(:)];
  end

  tol1(k) = abs(int_old(k)-int(k)); %absolute tolerance
  k       = find(tol1 > abs(tol*int)); %| tol1 > abs(tol));%indices to integrals which did not converge
   
  if any(k),% compute integrals again
      nk=length(k);%# of integrals we have to compute again
  else
    converge='y';
    break;
  end
  int_old=int;
end

int=reshape(int,N,M); % make sure int is the same size as the integration  limits
if nargout>1,
  tol1=reshape(tol1,N,M);
end

if converge=='n'
  if nk>1
    if (nk==N*M),
      disp('All integrals did not converge--singularities likely')
    else
      disp(sprintf('%d integrals did not converge--singularities likely', ...
		   nk))
    end
  else
    disp('Integral did not converge--singularity likely')
  end
end

if trace==1,
  plot(x_trace,y_trace,'+')
end


