% Module TRGAUSS in WAFO Toolbox. 
% Version 2.1.1   07-Sep-2005
% 
% 
% Readme        - Readme file for module TRGAUSS in WAFO Toolbox. 
%
% Misc
% createpdf     - PDF class constructor 
% pdfplot       - Plot contents of pdf structures 
% trplot        - Plot transformation, g, eg.  estimated with dat2tr.
%
%Transforms and non-linearities
% dat2gaus      - Transforms  x  using the transformation  g. 
% gaus2dat      - Transforms  xx  using the inverse of transformation  g. 
% mctrtest      - Test if stochastic process is Gaussian.
% spec2skew     - Estimate moments of 2'nd order waves due to Marthinsen and Winterstein 
% trangood      - Make transformation suitable for efficient transforms. 
% tranproc      - Transforms process X and up to four derivatives 
% trmak         - Put together a transformation object. 
% troptset      - Create or alter TRANSFORM OPTIONS structure. 
% trunmak       - Split a transformation object into its pieces.
%
% Transformed Gaussian model estimation
% cdf2tr        - Estimate transformation, g, from observed CDF. 
% dat2tr        - Estimate transformation, g, from data. 
% hermitetr     - Calculate transformation, g, proposed by Winterstein 
% ochitr        - Calculate transformation, g, proposed by Ochi et al. 
% lc2tr         - Estimate transformation, g, from observed crossing intensity. 
% lc2tr2        - Estimate transformation, g, from observed crossing intensity, version 2. 
%
%Gaussian probabilities and expectations
% bvnormcdf     - Bivariate Normal cumulative distribution function
% bvnormprb     - Bivariate Normal probability  
% mvnormpcprb   - Multivariate Normal probability with product correlation
% mvnormprb     - Multivariate Normal probability by Genz' algorithm.
% mvnortpcprb   - Multivariate Normal or student T probability with product correlation.
% rind          - Multivariate normal expectations
%
%Probability density functions (pdf) or intensity matrices
% chitwo2lc_sorm- SORM-approximation of crossing intensity for noncentral Chi^2 process
% chitwo2lc_sp  - Saddlepoint approximation of crossing intensity for noncentral Chi^2 process 
% dirsp2chitwo  - Parameters in non-central CHI-TWO process for directional Stokes waves.
% iter          - Calculates a Markov matrix fmM  given a rainflow matrix frfc; 
% iter_mc       - Calculates a kernel f_xy of a MC given a rainflow matrix 
% mc2rfc        - Calculates a rainflow matrix given a Markov chain with kernel f_xy; 
% mctp2rfc      - Calculates a rainflow matrix given a Markov matrix f_mM of a Markov chain 
% mctp2tc       - Calculates frequencies for the upcrossing troughs and crests 
% nt2fr         - Calculates the frequency matrix given the counting distribution matrix. 
% spec2cmat     - Joint intensity matrix for cycles 
% spec2mmtpdf   - Calculates joint density of Maximum, minimum and period. 
% spec2tccpdf   - Evaluates densities of wave period Tcc, wave lenght Lcc. 
% spec2thpdf    - Joint density of amplitude and period/wave-length characteristics 
% spec2tpdf     - Evaluates densities for crest-,trough-period, length. 
% spec2tpdf2    - Evaluates densities for various wave periods or wave lengths 
% specq2lc      - Saddlepoint approximation of crossing intensity for quadratic sea. 
% th2vhpdf      - Transform joint T-H density to V-H density 
%
%Cumulative distribution functions (cdf)
% lomaxcdf      - CDF for local maxima for a zero-mean Gaussian process 
% spec2acat     - Evaluates survival function R(h1,h2)=P(Ac>h1,At>h2). 
% spec2acdf     - Evaluates cdf of crests P(Ac<=h) or troughs P(At<=h). 
%
