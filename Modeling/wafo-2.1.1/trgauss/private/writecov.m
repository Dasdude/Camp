function writecov(cov,nr)
% WRITECOV Calculates spline coefficients for the covariance and 
%          its derivatives  r_X^(i)(t), i = 0,1,2,3,4.
%          The results are saved on ascii files  Cd*.in,
%          and are used by  minmax, wave_t  and  wave_th.
%
%  CALL: writecov(cov);
%        writecov(cov,nr);
% 
%        cov  = [T,r,r1,r2,...], as returned from  spec2cov.
%        nr   = number of derivatives to write
%               (default=all=size(cov,2)-2).

if nargin<2
  nr=size(cov,2)-2;
end

if size(cov,2)<2+nr
  error(['You must supply at least nr=', int2str(nr), ' derivatives.'])
end

for k=0:nr
  filename=['Cd', int2str(k), '.in'];
  if exist(filename)
    delete(filename)
  end
end

n=size(cov,1);
for k=0:nr
  filename=['Cd', int2str(k), '.in'];
  covar=[cov(:,1), cov(:,k+2), zeros(n,3)];
  fid=fopen(filename,'wt');
  fprintf(fid,'%12.10f %12.10E %4.2f %4.2f %4.2f\n',covar');
  fclose(fid);
end
