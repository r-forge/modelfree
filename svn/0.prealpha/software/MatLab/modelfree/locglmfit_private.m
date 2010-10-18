function [pfit,etafit,H] = locglmfit_private(xfit,r,m,h,x,link,guessing,...
    lapsing,K,p,ker,maxiter,tol)
%
% THIS IS AN INTERNAL FUNCTION: USE LOCGLMFIT FOR BEST RESULTS
%
% Fisher scoring method for local polynomial estimator of a psychometric
% function (PF). This function is a private function used by LOCGLMFIT. 
%
% INPUT
%
% xfit - points in which to calculate the estimate
% r - number of successes in points x
% m - number of trials in points x 
% h - bandwidths
% x - stimulus levels 
% link - name of the link function to be used
% guessing - guessing rate
% lapsing - lapsing rate
% K - power parameter for Weibull and reverse Weibull link
% p - degree of the polynomial 
% ker - kernel function for weights 
% maxiter - maximum number of iterations in Fisher scoring
% tol - tolerance level at which to stop Fisher scoring
%
% OUTPUT
%
% pfit - value of the local polynomial estimate in points x
% etafit - estimate of eta (link of pfit)
% H - Hat matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PROGRAM

%%%% INITAILS VALUES

% link
lims = [ guessing, 1 - lapsing ];
if ( strcmp( link, 'weibull' ) || strcmp( link, 'revweibull' ) )
     LINK = feval( strcat( link, '_link' ), K, lims );
else
     LINK = feval( strcat( link, '_link' ), lims );
end
linkfun = LINK{1};
derinv  = LINK{2};
linkinv = LINK{3};

n = length(r);
nx = length(xfit);

% data into column vectors
tmp1(:,1) = r;
r = tmp1;

tmp1(:,1) = m;
m = tmp1;

tmp1(:,1) = x;
x = tmp1;

tmp2(1,:) = xfit;
xfit = tmp2;

%%%%
%%%%  MATRICES FOR CALCULATION

% vector of differences xfit-x
diffx = reshape( repmat(xfit,n,1)-repmat(x,1,nx), n*nx,1);

% calculate weights given by the kernel ker and bandwidth h
kerx = zeros( nx * length( x ), 1 );
if length( h ) == 1
    kerx = feval(ker,diffx,0,h);
else
    h_mat = reshape( repmat(h',n,1), n*nx,1);
    kerx = feval(ker,diffx,0,h_mat);
end

% form a matrix of 1,x,x^2,...,x^p
tmpX0 = repmat(diffx,1,p+1).^repmat(0:p,n*nx,1);

% intial values for the algorithm
mu0 = (r + .5)./(m + 1);

eta0 = linkfun(mu0);
Z0 = repmat(eta0,nx,1);
eta_mu = feval(derinv,mu0);

% binomial weights for the algorithm
W = diag(repmat(m ./(eta_mu.*sqrt(mu0 .* (1 - mu0))),nx,1));

% kernel weights
KerX = diag(sqrt(kerx(:,1)));

% combined binomial and kernel weights
WK = W .* KerX;

% raw means in a matrix form
KM = repmat( r./m , nx , 1 );

% create X0 matrix
X0 = zeros(n*nx,(p+1)*nx);

for l = 1:nx,
    X0((1:n)+(l-1)*n,(1:p+1)+(l-1)*(p+1)) = tmpX0((1:n)+(l-1)*n,:);
end

% linear estimator
X = WK * X0;
Y = WK * Z0; 
beta = inv(X'*X)*(X'*Y);

% inital values for stopping the loop
iternum = 0;
etadiff = tol+1;
etafit = (X0 * beta);
mu_raw = repmat(mu0,nx,1);
M = repmat(m,nx,1);
score = 1;

% offset value (ensure no limiting values appear in the algorithm)
epsilon = 1e-5;

%%%%%
%%%%% FISHER SCORING
while ((iternum<maxiter)&&(etadiff>tol)&&(score)),

    % obtain values from previous loop
    mu_old = mu_raw;
    eta_old = etafit;
    
    % new mean
    mu = linkinv(etafit);
    mu_raw = mu;
    
    % offset
    mu(mu>=1-lapsing-epsilon) = 1-lapsing-epsilon;
    mu(mu<=guessing+epsilon) = guessing+epsilon;
    
    % derivatived d eta / d mu
    eta_mu = derinv(mu);

    % z scores
    z = etafit + (KM - mu) .* eta_mu;
    
    % new weights
    WK = diag(1./(eta_mu.*sqrt(mu .* (1 - mu)./M))) .* KerX;

    % linear estimator
    X = WK * X0;
    Y = WK * z;
    
    % new estimate of beta
    beta = (X'*X)\(X'*Y); %    beta = inv(X'*X)*(X'*Y);
    
    % beta0 (i.e. value of eta function)
    eta1 = beta(1:(p+1):end);

    % new estiate of eta and its derivatives
    etafit = (X0 * beta);
    
    % increase iteration count and adjust stopping values
    iternum = iternum + 1;
    mudiff = max(max(abs(mu_old-mu_raw)));
    etadiff = max(max(abs(eta_old-etafit)));
    score = ~((mudiff<tol)&(max(abs(eta1))>50));
end

% warning about exceeding iteration max
if (maxiter==iternum),
    warning('iteration limit reached')
end

% Hat matrix
if nargout == 3,
    tmpH = (X'* X)\(X' * WK);
    tmpH = tmpH(1:(p+1):end,:);
    H = zeros(nx,n);
    for i = 1:nx,
        H(i,:) = tmpH(i,(1:n)+(i-1)*n);
    end
end

% retrive beta0 and remove v. large and v. small values
etafit = beta(1:(p+1):end);

% find estimate of PF
pfit = feval(linkinv,etafit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KERNELS

% Epanechnikov
function y = epanechnikov(x,m,s)
X = x/s;
X(abs(X)>1) = 1;
y = 0.75 * (1 - X.^2)/s;

% triangular
function y = triangular(x,m,s)
X = x/s;
X(abs(X)>1) = 1;
y = (1 - abs(X))/s;

% tri-cube
function y = tricube(x,m,s)
X = x/s;
X(abs(X)>1) = 1;
y = ((1 - abs(X).^3).^3)/s;

% bi-square
function y = bisquare(x,m,s)
X = x/s;
X(abs(X)>1) = 1;
y = ((1 - abs(X).^2).^2)/s; 

% uniform
function y = uniform(x,m,s)
X = x/s;
y = (abs(X)<=1)/2;