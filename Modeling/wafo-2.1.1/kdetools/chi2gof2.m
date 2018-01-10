function [pvalue, test, v] = chi2gof2(fs,ft)
% CHI2GOF2 CHI Squared Goodness Of Fit test.
%
% CALL:  P = chi2gof2(fs, ft)
%
%   P      = P - value  
%   fs     = fitted PDF matrix
%            evaluated at observed points.
%   ft     = theoretical PDF matrix
%            evaluated at simulated points 
%
%  Large P value -> good fit
%  Small P value -> lesser fit
%
%  Example : Check how well rayleigh data can be described by N(0,1)
%   xs  = wraylrnd(1,500,1);
%   xs2 = wnormrnd(0,1,100000,1); 
%   p   = chi2gof2(wnormpdf(xs),wnormpdf(xs2));
%           
% See also  wchi2cdf, qlevels2


%Tested on: matlab 5.3
% History:
% revised pab feb2005
% -updated see also line  
% revised pab 11.11.2000
% - made it independent of stats toolbox only dependent on wstats
% by pab 21.09.99  
  


% CHI^2 goodness of fit (GOF)   test

 if isstruct(fs) % structure
     r2=fs.f(:);
   else
    r2 = fs(:);
  end
ntresh=length(r2);
if ntresh>120 % only valid for samples larger than 120
   if isstruct(ft)
     r=ft.f;
   else
    r = ft;
  end
  k=max([ceil(sqrt(ntresh)),8]);%divide the data into k subsets (greater than 8)
  pk=100/k;                     % with equal probabilty
  np=ntresh/k; %the expected number of points in each group (must be greater than 5)
  
  grpdiv=qlevels2(r,(100-pk):-pk:pk); % find dividing levels 
  %grpdiv=normpdf(norminv( (pk:pk:(100-pk))/200))' 
 

  test=(sum(grpdiv(1 )>=r2)-np)^2/np +(sum(grpdiv(k-1)<r2)-np)^2/np; 
  for ix=2:k-1,
    test=test+(sum((grpdiv(ix-1 )<r2).*(grpdiv(ix )>=r2))-np)^2/np;
  end
  test %test statistic
  pvalue=1-wchi2cdf(test,k-1); % pvalue
  v=k-1;
else  
  pvalue=[]
  disp('to few data')
  ntresh
end
return


