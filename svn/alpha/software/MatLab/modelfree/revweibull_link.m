function link = revweibull_link(K,lims)
%
% Reverse Weibull link in cell form for use with GLMFIT, GLMVAL and other Matlab
% GLM functions 
%
% The guessingrate and 1-lapsing rate are fixed to values given in lims, and
% exponent is set to be equal K, hence link is a function of only one
% variable.   
%
% INPUT
% 
% K - exponent used in the definition of reverse Weibull link function
%
% OPTIONAL INPUT
%
% lims - two column vector specifying guessing rate and 1-lapsing rate; default
% is [0,1] 
%
% OUTPUT
%
% link - reverse Weibull link for use in all GLM functions; cell with 3 entries:  
%   	revweibullkFL - link function
%       revweibullkFD - derivative   
%   	revweibullkFI - inverse link

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET LINK
link = cell(3,1);
link{1} = @(x) revweibullkFL(x,g,l,K);
link{2} = @(x) revweibullkFD(x,g,l,K);
link{3} = @(x) revweibullkFI(x,g,l,K);

% % % % % % % % % % % % % % % % % % 
% % % INTERNAL FUNCTIONS % % % % % 
% % % % % % % % % % % % % % % % % % 

%%%%%%%%%%%%%%%%% LINK FUNCTION DEFINITIONS%%%%%%%%%%%%%%%%%%%%%%
% REVERSE WEIBULLK
function eta = revweibullkFL(mu,g,l,k)
mu = max(min(l-eps,mu),g+eps);
eta = -(-log((mu-g)/(l-g))).^(1/k);

function eta = revweibullkFD(mu,g,l,k)
mu = max(min(l-eps,mu),g+eps);
eta = 1./(k*(-log((mu-g)/(l-g))).^((k-1)/k).*(mu-g));

function mu = revweibullkFI(eta,g,l,k)
eta = max(-(-log(eps/(l-g))).^(1/k),eta);
eta = min(-(-log(1-eps/(l-g))).^(1/k),eta);
mu = g + (l-g) * exp(-(-eta).^k);