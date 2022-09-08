# % Construct spherical grid over (phi,theta) based on input parameters:
# % - nOctants: number of "octants"; for each increment in theta, nOctants
# %   additional points are added along phi; special cases: nOctants=0 and
# %   nOctants=-1
# % - phimax: largest value of phi (radians)
# % - GridSize: number of orientations between theta=0 and theta=pi/2
# % - closedPhi: set to true if grid point at phimax should be included
function powderpoints(exp)
  if exp.nOctants == 1
    phimax = pi*0.25
  elseif exp.nOctants == 2
    phimax = pi*0.5
  elseif exp.nOctants == 3
    phimax = pi*0.75
  elseif exp.nOctants == 4
    phimax = pi
  end

  dtheta = 0.5*pi/exp.nResolution
  theta(j) = 0.5*pi*(1-(j-1)/exp.nResolution)
  dphi(j) = 0.5*pi/(j-1)
  phi(i,j) = 0.5*pi*(i-1)/(j-1)

  nPoints = 4*exp.nResolution^2 + 2

  Phi = zeros(nPoints)
  Phi[1] = 0.
  Theta = zeros(nPoints)
  Theta[1] = 0.
  Weight = zeros(nPoints)
  Weight[1] = 0.

  
  0.5*pi*(i-1)/(j-1)

  
  nOrientations = exp.nPoints + exp.nOctants*exp.nPoints*(exp.nPoints-1)*0.5
  phi = zeros(1,nOrientations)
  theta = zeros(1,nOrientations)
  Weights = zeros(1,nOrientations)
  sindth2 = sind(dtheta*0.5)
  w1 = 0.5

  #North pole (z orientation)
  phi[1] = 0
  theta[1] = 0
  Weights[1] = phimax*(1-cosd(dtheta*0.5))

  #All but equatorial slice
  phistart = 2
  for ii = 2:exp.nPoints-1
    nphi = nOct*(ii-1)+1
    dPhi = phimax/(nPhi-1)
    idx = Start+(0:nPhi-1)
    Weights[idx] = 2*sin((ii-1)*dtheta)*sindth2*dPhi*[w1 ones(1,nPhi-2) .5]
    phi[idx] = linspace(0,phimax,nPhi)
    theta[idx] = (ii-1)*dtheta
    phistart = phistart + nPhi
  end

  #   % Equatorial slice
  nPhi = nOct*(exp.nPoints-1)+1;
  dPhi = phimax/(nPhi-1);
  idx = phistart + (0:nPhi-1);
  phi(idx) = linspace(0,phimax,nPhi);
  theta(idx) = pi/2;
  Weights(idx) = sindth2*dPhi*[w1 ones(1,nPhi-2) 0.5];

# #   % Border removal
#   if ~closedPhi
#     rmv = cumsum(nOct*(1:exp.nPoints-1)+1)+1;
#     phi(rmv) = [];
#     theta(rmv) = [];
#     Weights(rmv) = [];
#   end

Weights = 2*(2*pi/phimax)*Weights; #% sum = 4*pi

end
end