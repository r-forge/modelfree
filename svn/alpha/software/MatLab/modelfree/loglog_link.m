function link = loglog_link( lims )
%
% Loglog link in cell form for use with GLMFIT, GLMVAL and other Matlab
% GLM functions 
%
% The guessing rate and 1-lapsing rate are fixed to values given in lims
% hence link is a function of only one  variable.    
%
% OPTIONAL INPUT
%
% lims - two column vector specifying guessing rate and 1-lapsing rate; default
% is [0,1] 
%
% OUTPUT
%
% link - Loglog link for use in all GLM functions; cell with 3 entries:  
%   	loglogFL - link function
%       loglogFD - derivative   
%   	loglogFI - inverse link
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
link{1} = @(x) loglogFL( x, g, l );
link{2} = @(x) loglogFD( x, g, l );
link{3} = @(x) loglogFI( x, g, l );

% % % % % % % % % % % % % % % % % % 
% % % INTERNAL FUNCTIONS % % % % % 
% % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%
% LOGLOG WITH LIMITS

% link
function eta = loglogFL(mu,g,l)

mu = max(min(l-eps,mu),g+eps);
eta = -log(-log((mu-g)/(l-g)));

% derivative
function eta = loglogFD(mu,g,l)

mu = max(min(l-eps,mu),g+eps);
eta = -1./((mu-g).*log((mu-g)./(l-g)));

% inverse link
function mu = loglogFI(eta,g,l)

eta = max(-log(-log((eps)/(l-g))),eta);
eta = min(-log(-log(1-eps/(l-g))),eta);
mu = g + (l-g) * exp(-exp(-eta));