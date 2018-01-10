function [int, tol1,k] = gaussq2d(fun,xlow,xhigh,ylow,yhigh,tol,p1,p2,p3,p4,p5,p6,p7,p8,p9)
%GAUSSQ2D Numerically evaluates a 2D integral using Gauss quadrature.
%
% CALL:  [int tol] = gaussq2d('Fun',xlow,xhigh,ylow,yhigh)
%      or
%        [int tol] = gaussq2d('Fun',xlow,xhigh,ylow,yhigh,reltol,p1,p2,...)
%
%       int = evaluated integral
%       tol = absolute tolerance abs(int-intold)
%       Fun = Fun(x,y), string containing the name of the function or  
%             a directly given expression enclosed in parenthesis.
%             The function may depend on parameters pi.
%xlow,xhigh = integration limits for x 
%ylow,yhigh = integration limits for y (obs. integration limits can
%                                             be vectors.)
%    reltol = relative tolerance (default 1e-3)
%
% Note that if there are discontinuities the region of integration 
% should be broken up into separate pieces.  And if there are singularities,
% a more appropriate integration quadrature should be used 
% (such as the Gauss-Chebyshev for a specific type of singularity).
%
%Example:%  integration for (x,y) in [0,1]x[-1,3] and [2,4]x[1,2]:
%
%   p1=2.; p2=0.5;
%   gaussq2d('(p2*x.^2.*y+p1)',[0 2],[1 4],[-1 1],[3 2],[],p1,p2)
%
% See also  gaussq, qrule2d

% tested on: matlab 5.2
% history:
% revised pab Nov2004
% -wrong input to comnsize fixed  
% -enabled function handle as input  
% revised pab 12Nov2003
% replaced call to distchk with comnsize
%modified by Per A. Brodtkorb 17.11.98 : 
% -accept multiple integrations limits
% -optimized by saving the weights as global constants and only  
%   computing the integrals which did not converge 
% -enabled the integration of directly given functions enclosed in
%   parenthesis 
% -adopted from NIT toolbox changed name from quad2dg to gaussq2d  

global bpx2 bpy2 wfxy2;
fv = [];
if isempty(bpx2)| isempty(bpy2)| isempty(wfxy2)
  [bpx2,bpy2,wfxy2] = qrule2d(2,2);
end

if exist('tol')~=1,
  tol=1e-3;
elseif isempty(tol),
  tol=1e-3;
end
[errorcode ,xlow,xhigh,ylow,yhigh]=comnsize(xlow,xhigh,ylow,yhigh);
if errorcode > 0
  error('Requires non-scalar arguments to match in size.');
end

[N M]=size(xlow);
%setup mapping parameters &  make sure they all are column vectors
xlow=xlow(:); xhigh=xhigh(:);ylow=ylow(:);yhigh=yhigh(:); 
nk=N*M;

%Map to x
qx=(xhigh-xlow)/2;
px=(xhigh+xlow)/2;
x=permute(qx(:,ones(1,2), ones(1,2) ),[ 2 3 1]).*bpx2(:,:,ones(1,nk)) ...
    +permute(px(:,ones(1,2), ones(1,2) ),[2 3 1]);
%Map to y
qy=(yhigh-ylow)/2;
py=(yhigh+ylow)/2;
y=permute(qy(:,ones(1,2), ones(1,2) ),[ 2 3 1]).*bpy2(:,:,ones(1,nk)) ...
    +permute(py(:,ones(1,2), ones(1,2) ),[2 3 1]);

if (isa(fun,'char') & any(fun=='(')),
  exec_string=['fv=', fun, ';']; %the call function is already setup
else 
  %set up the call function  
  if isa(fun,'function_handle')
    fun = func2str(fun);
  end
  exec_string=['fv=', fun, ' (x,y'];
  num_parameters=nargin-6;
  for i=1:num_parameters,
    exec_string=[exec_string,',p',int2str(i)];
  end
  exec_string=[exec_string,');'];
end
eval(exec_string);
int_old = squeeze(sum(sum(wfxy2(:,:,ones(1,nk)).*fv))).*qx.*qy;
int=zeros(size(int_old));tol1=int;
k=1:nk;

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to int)
%  2) the last update is very small (compared to tol)
%  3) There are more than 8 iterations. This should NEVER happen. 

converge='n';
for i=1:8,
  gnum=int2str(2^(i+1));
  eval(['global bpx',gnum,' wfxy',gnum,'  bpy',gnum,';']);
  if isempty(eval(['bpx',gnum]))| isempty(eval(['bpy',gnum])) ,  
    eval(['[bpx',gnum,',bpy',gnum,',wfxy',gnum,']=qrule2d(',gnum,',', gnum,');']);
  end
  eval(['x=permute(qx(k,ones(1,',gnum,'), ones(1,',gnum,') ),[ 2 3 1]).*bpx',gnum,'(:,:,ones(1,nk)) +permute(px(k,ones(1,',gnum,'), ones(1,',gnum,') ),[2 3 1]);']);

  eval(['y=permute(qy(k,ones(1,',gnum,'), ones(1,',gnum,') ),[ 2 3 1]).*bpy',gnum,'(:,:,ones(1,nk)) +permute(py(k,ones(1,',gnum,'), ones(1,',gnum,') ),[2 3 1]);']);
    eval(exec_string);
    eval(['int(k) = squeeze(sum(sum((wfxy',gnum,'(:,:,ones(1,nk)).*fv)))).*qx(k).*qy(k);']);
    tol1(k)=abs(int_old(k)-int(k)); % absolute tolerance
    k=find(tol1 > abs(tol*int)|tol1 > abs(tol));%indices to integrals which did not converge
  
    if any(k), % compute integrals again
      nk=length(k);%# of integrals we have to compute again
    else
      converge='y';
      break;
    end
    
     int_old(k)=int(k);
end
int=reshape(int,[N M]); % reshape int to the same size as input arguments
if converge=='n',
  disp('Integral did not converge--singularity likely')
end



