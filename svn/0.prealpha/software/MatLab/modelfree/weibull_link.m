function link = weibull_link(K,lims)
%
% Weibull link in cell form for use with GLMFIT, GLMVAL and other Matlab
% GLM functions 
%
% The guessing rate and 1-lapsing rate are fixed to values given in lims, and
% exponent is set to be equal K, hence link is a function of only one
% variable.    
%
% INPUT
% 
% K - exponent used in the definition of Weibull link function
%
% OPTIONAL INPUT
%
% lims - two column vector specifying guessing rate and 1-lapsing rate; default
% is [0,1] 
%
% OUTPUT
%
% link - Weibull link for use in all GLM functions; cell with 3 entries:  
%   	weibullkFL - link function
%       weibullkFD - derivative   
%   	weibullkFI - inverse link

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PROGRAM

%%%%
%%%% CHECK INPUT PARAMETERS + INFORM OF DEFAULT VALUES
if (nargin<1)
    error('Exponent of reverse Weibull is mandatory');
end
if (nargin<2), 
    g=0; l=1;
else 
    g = min(lims);
    l = max(lims);
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
checkinput( 'exponentk', K );
checkinput( 'guessingandlapsing', [ g, l ] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET LINK
link = cell(3,1);
link{1} = @(x) weibullkFL(x,g,l,K);
link{2} = @(x) weibullkFD(x,g,l,K);
link{3} = @(x) weibullkFI(x,g,l,K);

% % % % % % % % % % % % % % % % % % 
% % % INTERNAL FUNCTIONS % % % % % 
% % % % % % % % % % % % % % % % % % 

%%%%%%%%%%%%%%%%% LINK FUNCTION DEFINITIONS%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% WEIBULLK
function eta = weibullkFL(mu,g,l,k)
mu = max(min(l-eps,mu),g+eps);
eta = (-log((l-mu)./(l-g))).^(1/k);

function eta = weibullkFD(mu,g,l,k)
mu = max(min(l-eps,mu),g+eps);
eta = 1./(k*(-log((l-mu)./(l-g))).^((k-1)/k).*(l-mu));

function mu = weibullkFI(eta,g,l,k)
eta = max((-log(1-eps./(l-g))).^(1/k),eta);
eta = min((-log(eps./(l-g))).^(1/k),eta);
mu = g + (l-g) * (1 - exp(-eta.^k));