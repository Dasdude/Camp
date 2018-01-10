function Snew=rotspec(S,phi,rotateGrid)
%ROTSPEC Rotate spectrum anti-clockwise around the origin. 
%
% CALL:  Snew=rotspec(S,phi,rotateGrid)
%
%       Snew = new rotated spectrum
%          S = spectrum
%        phi = rotation angle (default 0)
% rotateGrid = 1 if rotate grid of Snew physically (thus Snew.phi=0).
%              0 if rotate so that only Snew.phi is changed  
%                (the grid is not physically rotated)  (default)
%
% ROTSPEC Rotates a spectrum anti-clockwise around the origin.
% The spectrum can be of any of the two-dimensional types.
% For spectrum in polar representation:
%    newtheta = theta+phi, but circulant such that -pi<newtheta<pi
% For spectrum in Cartesian representation:
%    If the grid is rotated physically, the size of it is preserved
%    (maybe it must be increased such that no nonzero points are
%    affected, but this is not implemented yet: i.e. corners are cut off)
% The spectrum is assumed to be zero outside original grid.
% NB! The routine does not change the type of spectrum, use spec2spec
%     for this.
%
% Example
%  S=demospec('dir');
%  wspecplot(S), hold on  
%  wspecplot(rotspec(S,pi/2),'r'), hold off
%  
% See also datastructures, spec2spec

% History:   
% Revised pab mar 2005
% Revised pab Jan 2005
% fixed an illegal else statement
% Revised pab Sept 2004:
%  -Removed old unused code
%  -Added: Example and rotateGrid to input  
% revised IR Aug 2004  
% revised IR, PAB
%  - BUG: Made sure  -pi<newtheta<pi
% by es 17.07.1999
  
  
% TODO % Make physical grid rotation of cartesian coordinates more robust.

error(nargchk(1,3,nargin))
if (nargin<2 | isempty(phi))
  phi =0;
end
if (nargin<3 | isempty(rotateGrid))
  rotateGrid = 0;
end

Snew=S;

if (~isfield(S,'phi') | isempty(S.phi)),
  Snew.phi=0; 
end
Snew.phi=mod(Snew.phi+phi+pi,2*pi)-pi;

switch lower(S.type(end-2:end))
  case 'dir',
   % any of the directinal types 
   % Make sure theta is from -pi to pi
   phi0       = Snew.theta(1)+pi; 
   Snew.theta = Snew.theta-phi0;
   
   % make sure -pi<phi<pi
   Snew.phi   = mod(Snew.phi-phi0+pi,2*pi)-pi;
   if (rotateGrid & (Snew.phi~=0))
     % Do a physical rotation of spectrum
     theta = Snew.theta;
     ntOld = length(theta)
     nt = ntOld-(theta(1)==theta(end));
     Snew.theta(1:nt) = mod(Snew.theta(1:nt)+Snew.phi+pi,2*pi)-pi;
     Snew.phi = 0;
     [Snew.theta,ind] = sort(Snew.theta(:));
     Snew.S = Snew.S(ind,:);
     
     if (Snew.theta(1)==-pi)
        if (nt<ntOld)
          Snew.S(ntOld,:) = Snew.S(1,:);
        else
           Snew.S(nt+1,:) = Snew.S(1,:);
        end
     elseif (nt<ntOld)
        Snew.S(ntOld,:)   = []
        Snew.theta(ntOld) = [];
     end
   end
    
 case 'k2d', 
  % any of the 2D wave number types
   %Snew.phi   = mod(Snew.phi+phi+pi,2*pi)-pi;
  
  
  if (rotateGrid & (Snew.phi~=0))
    % Do a physical rotation of spectrum
    
    [k,k2] = meshgrid(S.k,S.k2);
    [th,r] = cart2pol(k,k2);
    [k,k2] = pol2cart(th+Snew.phi,r);
    Sn = interp2(S.k,S.k2,S.S,k,k2);
    Sn(isnan(Sn))=0.;
    Snew.S   = Sn;
    Snew.phi = 0;
  end
otherwise
  %disp('Can only rotate two dimensional spectra')
end
  
return

