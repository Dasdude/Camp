function [f_mM, paramu] = wminmax(spec,nit,paramu,t)
%WMINMAX Calculates joint density of minimum and following maximum
%        in a zero-mean stationary Gaussian with normalized spectrum spec
%
% CALL:  f_mM   = wminmax(spec,nit,paramu,t);     
%
%        f_mM  = joint density of minimum and following maximum in  X(t),
%                for def=0 otherwise f_hh is the joint density of (H_1,H_2)
%                that is heights of crest and trough in up-crossing waves.
%                
%        ALL THIS INPUTS NEEDS TO BE SET (NO DEFAULT VALUES ARE ALLOWED)
%    spec   = normalized spectrum L0=L2=1
%    nit    = order of numerical integration: 0,1,2,3,4,5.
%    paramu = parameter vector defining discretization of min/max values.
%    t      = grid of time points between maximum and minimum (to
%             integrate out). interval between maximum and the following
%             minimum,
  
%History
% revised pab Dec2003
% -  replaced code with call to spec2cov2
paramv = paramu;
par    = [paramu, paramv];


IAC = 1;
ISQ = 0;
EPS = 5e-5;
EPSS = 1e-6;
EPS0 = 1e-5;
tol = [IAC,ISQ,EPS,EPSS,EPS0];
%tol=[1, 0, 5e-5, 1e-6, 1e-5];

%g=[(-5:0.02:5)', (-5:0.02:5)']; 
 
g = spec.tr;

if length(t)>101
  error('nr. of time points limited to 101.')
end


if abs(t(1))>0.00001 
  error('t(1) < or > 0.')
end
if length(t) < 2
  error('nr. of wavelength <2.')
end

if par(3)<1
  error('Require n>0.')
end
if par(6)<1
  error('Require n>0.')
end

accur = [nit tol];

dt    = t(2)-t(1);
if 1,
  Nt = length(t)-1;
  nr = 4;
  R = spec2cov2(spec,nr,Nt,dt);
  cov = [t(:) R];
else
  S1 =specinterp(spec,dt);
  R = spec2cov(S1,4,400,4);
  cov = [R.t, R.R, R.Rt, R.Rtt,  R.Rttt, R.Rtttt];
end

writecov(cov,4);

if exist('t.in'), delete('t.in'), end
if exist('transf.in'), delete('transf.in'), end
if exist('accur.in'), delete('accur.in'), end
if exist('Mm.in'), delete('Mm.in'), end

disp('   Writing data.')
fid = fopen('t.in','wt');
fprintf(fid,'%8.5f\n',t);
fclose(fid);

fid = fopen('accur.in','wt');
fprintf(fid,'%2.0f %2.0f %2.0f\n',accur(1:3));
fprintf(fid,'%8.7e %8.7e %8.7e\n',accur(4:6));
fclose(fid);

fid = fopen('transf.in','wt');
fprintf(fid,'%8.5e %8.5e \n',g');
fclose(fid);

fid=fopen('Mm.in','wt');
fprintf(fid,'%8.6e %8.6e %3.0f\n',par(1:3));
fprintf(fid,'%8.6e %8.6e %3.0f\n',par(4:6));
fclose(fid);

disp('   Starting Fortran executable.')
dos([wafoexepath, 'minmax.exe']);

disp('   Loading data.')
f_mM = loaddata('Maxmin.out');
f_mM = reshape(f_mM(:,3),paramv(3),paramu(3));
f_mM = rot90(f_mM,-2);





