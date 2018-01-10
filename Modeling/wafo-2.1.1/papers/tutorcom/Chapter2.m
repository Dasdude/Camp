%% CHAPTER2 Modelling random loads and stochastic waves
%
% Chapter2 contains the commands used in Chapter 2 of the tutorial and
% present some tools for analysis of randodm functions with
% respect to their correlation, spectral and distributional properties.
% The code is divided into three examples: 
%
% Example1 is devoted to estimation of different parameters in the model.
% Example2 deals with spectral densities and
% Example3 presents the use of WAFO to simulate samples of a Gaussian
% process.
%
% Some of the commands are edited for fast computation. 
% Each set of commands is followed by a 'pause' command.
% 

% Tested on Matlab 5.3, 7.01
% History
% Revised pab sept2005
%  Added sections -> easier to evaluate using cellmode evaluation.
% Revised pab Dec2004
% Created by GL July 13, 2000
% from commands used in Chapter 2
%

pstate =  'off';

%% Section 2.1 Introduction and preliminary analysis
%% Example 1: Sea data
%% Observed crossings compared to the expected for Gaussian signals
xx = load('sea.dat');
me = mean(xx(:,2))
sa = std(xx(:,2))
xx(:,2) = xx(:,2) - me;
lc = dat2lc(xx);
pause(pstate)


%% Turningpoints and irregularity factor
T = max(xx(:,1))-min(xx(:,1))
f0 = interp1(lc(:,1),lc(:,2),0)/T
pause(pstate)

tp = dat2tp(xx);
alfa = f0/(length(tp)/(2*T))
pause(pstate)

% A part of sea data is visulized with the following commands
clf
waveplot(xx,tp,'k-','*',1,1)
axis([0 2 -inf inf])
wafostamp([],'(ER)')
pause(pstate)

%% Finding possible spurious points
%However, if the amount of data is too large for visual examinations one
%could use the following criteria to find possible spurious points. One
%must be careful using the criteria for extremevalue analysis, because
%these criteria might remove the highest and steepest waves.
dt = diff(xx(1:2,1));
dcrit = 5*dt;
ddcrit = 9.81/2*dt*dt;
zcrit = 0;
[inds, indg] = findoutliers(xx,zcrit,dcrit,ddcrit);
pause(pstate)

%% Section 2.2 Frequency Modeling of Load Histories
%% Periodogram: Raw spectrum
clf
S = dat2spec2(xx,9500);
wspecplot(S)
wafostamp([],'(ER)')
pause(pstate)

%% Calculate moments  
mom = spec2mom(S,4)
[sa sqrt(mom(1))]
pause(pstate)

%% Section 2.2.1 Random functions in Spectral Domain - Gaussian processes
%% Smoothing of spectral estimate 
clf
S1 = dat2spec2(xx,200);
S2 = dat2spec2(xx,50);
wspecplot(S1,[],'-.')
hold on
wspecplot(S2)
hold off
wafostamp([],'(ER)')
pause(pstate)

%% Estimated autocovariance
clf
R2 = spec2cov(S1,1);
Rest = dat2cov(xx,80,[],'- -');
covplot(R2,80,[],'.')
hold on
covplot(Rest)
wafostamp([],'(ER)')
hold off
pause(pstate)

%% Section 2.2.2 Transformed Gaussian models
rho3 = wskewness(xx(:,2))
rho4 = wkurtosis(xx(:,2))

[sk, ku]=spec2skew(S1)

%% Comparisons of 3 transformations
clf
gh = hermitetr([],[sa sk ku me]);
g  = gh; g(:,2)=g(:,1)/sa;
trplot(g)

[glc, test0] = dat2tr(xx);
hold on
plot(glc(:,1),glc(:,2),'b-')
plot(gh(:,1),gh(:,2),'b-.')
hold off
wafostamp([],'(ER)')
pause(pstate)

%%  Test Gaussianity of a stochastic process.
%MCTRTEST simulates  e(g(u)-u) = int (g(u)-u)^2 du  for Gaussian processes 
%  given the spectral density, S. The result is plotted if test0 is given.
%  This is useful for testing if the process X(t) is Gaussian.
%  If 95% of TEST1 is less than TEST0 then X(t) is not Gaussian at a 5% level.
      
%the following test takes time
N = length(xx);
test1 = mctrtest(S1,[N,50],test0);
wafostamp([],'(CR)')
pause(pstate)

%% Normalplot of data xx
% indicates that the underlying distribution has a "heavy" upper tail and a
% "light" lower tail. 
clf
wnormplot(xx(:,2))
wafostamp([],'(ER)')
pause(pstate)

%% Section 2.2.3 Spectral densities of sea data
%% Example 2: Different forms of spectra
clf
Hm0 = 7; Tp = 11;
spec = jonswap([],[Hm0 Tp]);
spec.note
pause(pstate)


%% Directional spectrum and Encountered directional spectrum
clf
D = spreading(101,'cos2s',0,[],spec.w,1)
Sd = mkdspec(spec,D)
pause(pstate)

clf
Se = spec2spec(Sd,'encdir',0,10);
wspecplot(Se), hold on
wspecplot(Sd,1,'--'), hold off
wafostamp([],'(ER)')
pause(pstate)

%% Frequency spectra
clf
S1 =spec2spec(Sd,'freq');
S2 = spec2spec(Se,'enc');
wspecplot(spec), hold on
wspecplot(S1,1,'.'),
wspecplot(S2),
wafostamp([],'(ER)')
hold off
pause(pstate)

%% Wave number spectrum
clf
Sk = spec2spec(spec,'k1d')
Skd = spec2spec(Sd,'k1d')
wspecplot(Sk), hold on
wspecplot(Skd,1,'--'), hold off
wafostamp([],'(ER)')
pause(pstate)

%% Effect of waterdepth on spectrum
clf
wspecplot(spec,1,'--'), hold on
S20 = spec;
S20.S = S20.S.*phi1(S20.w,20);
S20.h = 20;
wspecplot(S20),  hold off
wafostamp([],'(ER)')
pause(pstate)

%% Section 2.3 Simulation of transformed Gaussian process
%% Example 3: Simulation of random sea    
% Reconstruct replaces the spurious points of seasurface by simulated
% data on the basis of the remaining data and a transformed Gaussian
% process. As noted previously one must be careful using the criteria 
% for finding spurious points when reconstructing a dataset, because
% these criteria might remove the highest and steepest waves as we can see
% in this plot where the spurious points is indicated with a '+' sign:
%

clf
[y, grec] = reconstruct(xx,inds);
waveplot(y,'-',xx(inds,:),'+',1,1)
axis([0 inf -inf inf])
wafostamp([],'(ER)')
pause(pstate)

%%
clf
L = 200;
x = dat2gaus(y,grec);
Sx = dat2spec(x,L);
pause(pstate)
      
clf
dt = spec2dt(Sx)
Sx.tr = grec;
ysim = spec2sdat(Sx,480);
waveplot(ysim,'-')
wafostamp([],'(CR)')
pause(pstate)
 
%% Estimated spectrum compared to Torsethaugen spectrum
clf
Tp = 1.1;
H0 = 4*sqrt(spec2mom(S1,1))
St = torsethaugen([0:0.01:5],[H0  2*pi/Tp]);
wspecplot(S1)
hold on
wspecplot(St,[],'-.')
wafostamp([],'(ER)')
pause(pstate)

%%
clf
Snorm = St;
Snorm.S = Snorm.S/sa^2;
dt = spec2dt(Snorm)
pause(pstate)

clf
[Sk Su] = spec2skew(St);
sa = sqrt(spec2mom(St,1));
gh = hermitetr([],[sa sk ku me]);
Snorm.tr = gh;
pause(pstate)

%% Transformed Gaussian model compared to Gaussian model
clf
dt = 0.5;
ysim_t = spec2sdat(Snorm,240,dt);
xsim_t = dat2gaus(ysim_t,Snorm.tr);
pause(pstate)

clf
xsim_t(:,2) = sa*xsim_t(:,2);
waveplot(xsim_t,ysim_t,5,1,sa,4.5,'r.','b')
wafostamp([],'(CR)')
pause(pstate)

