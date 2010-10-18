function link = probit_link( lims )
%
% Probit link in cell form for use with GLMFIT, GLMVAL and other Matlab
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
% link - Probit link for use in all GLM functions; cell with 3 entries:  
%   	probitFL - link function
%       probitFD - derivative   
%   	probitFI - inverse link
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
link{1} = @(x) probitFL( x, g, l );
link{2} = @(x) probitFD( x, g, l );
link{3} = @(x) probitFI( x, g, l );

% % % % % % % % % % % % % % % % % % 
% % % INTERNAL FUNCTIONS % % % % % 
% % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%
% PROBIT WITH LIMITS

% link
function eta = probitFL(mu,g,l)

mu = max(min(l-eps,mu),g+eps);
eta = norminv((mu-g)./(l-g));

% derivative
function eta = probitFD(mu,g,l)

mu = max(min(l-eps,mu),g+eps);
eta = 1./normpdf(norminv((mu-g)/(l-g)))/(l-g);

% inverse link
function mu = probitFI(eta,g,l)

eta = max(norminv(eps./(l-g)),eta);
eta = min(-norminv(eps./(l-g)),eta);
mu = g + (l-g) * normcdf(eta);