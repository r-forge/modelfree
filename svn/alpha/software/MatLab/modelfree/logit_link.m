function link = logit_link( lims )
%
% Logit link in cell form for use with GLMFIT, GLMVAL and other Matlab
% GLM functions 
%
% The guessing rate and 1-lapsing rate are fixed to values given in lims hence 
% link is a function of only one variable.    
%
% OPTIONAL INPUT
%
% lims - two column vector specifying guessing rate and 1-lapsing rate; default
% is [0,1] 
%
% OUTPUT
%
% link - Logit link for use in all GLM functions; cell with 3 entries:  
%   	logitFL - link function
%       logitFD - derivative   
%   	logitFI - inverse link
%
% Created by Ivan Marin-Franch, 20/03/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PROGRAM
if (nargin<1), 
    g=0; l=1;
else 
    g = min(lims);
    l = max(lims);
end
checkinput( 'guessingandlapsing', [ g, l ] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET LINK
link = cell(3,1);
link{1} = @(x) logitFL( x, g, l );
link{2} = @(x) logitFD( x, g, l );
link{3} = @(x) logitFI( x, g, l );

% % % % % % % % % % % % % % % % % % 
% % % INTERNAL FUNCTIONS % % % % % 
% % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%
% LOGIT WITH LIMITS

% link
function eta = logitFL(mu,g,l)

mu = max(min(l-eps,mu),g+eps);
eta = log((mu-g)./(l-mu));

% derivative
function eta = logitFD(mu,g,l)

mu = max(min(l-eps,mu),g+eps);
eta = (l-g)./((mu-g).*(l-mu));

% inverse link
function mu = logitFI(eta,g,l)

eta = max(log(eps./(l-g)),eta);
eta = min(log((l-g)./eps),eta);
mu = g + (l-g)./(1+exp(-eta));