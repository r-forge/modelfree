function pfit = binomval_lims( b, xfit, link, lims, K)
%
% Fitted values at points xfit for a binomial GLM with coefficients b and
% lower and upper constrains as in lims.  
%
% INPUT
%
% xfit - points in which to calculate the estimate
% b - column vector of coefficients (result of BINOFIT_LIMS)
%
% OPTIONAL INPUT
%
% link - link function; choose from: 'logit', 'probit', 'loglog',
% 'comploglog', 'weibull', 'revweibull'; should be the same as used in
% BINOFIT_LIMS but default is the canonical link   
% lims - vector of length 2 with guessing and 1-lapsing rate; should be the
% same as used in BINOFIT_LIMS but default is [0,1]
% K - the exponent for weibull and revweibull link function; default is 2
%
% OUTPUT
%
% pfit - fitted values; column vector of the same length as xfit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PROGRAM

%%%% CHECK INPUT PARAMETERS + INFORM OF DEFAULT VALUES
% First 2 paramaters are mandatory
if (nargin<2)
    error('Coefficients (b) and stimulus levels (xfit) are mandatory');
end

%%%% DEFAULTS
if (nargin<3)
    link = 'logit';
    disp('default link function is ''logit''');
end

if (nargin<4)
    lims = [0 1];
    disp('default lower and upper limits are 0 and 1');
end

if (nargin<5)
    K = 2;
    if strcmp(link, 'weibull')
        disp('default exponent for Weibull link function is 2');
    elseif strcmp(link, 'revweibull')
        disp('default exponent for reverse Weibull link function is 2');
    end
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
checkinput( 'designpoints', xfit );
checkinput( 'linkfunction', link );
checkinput( 'guessingandlapsing', lims );
checkinput( 'exponentk', K );

%%%%
%%%% INITIALS VALUES

% column vectors
tmp1(:,1) = xfit;
xfit = tmp1;

% initial values
Lxfit = length(xfit);
p = length(b)-1;

% create matrix with all powers (1 to p) of xfit (for use in GLMFIT)
x = repmat(xfit,1,p).^repmat((1:p),Lxfit,1);

%%%%%
% assign link function
switch link
case'logit' 
% LOGISTIC
    linkfun = logit_link( lims );
case'probit' 
% PROBIT
    linkfun = probit_link( lims );
case'loglog'
% LOG-LOG
    linkfun = loglog_link( lims );
case'comploglog' 
% COMPLEMENTARY LOG-LOG
    linkfun = comploglog_link( lims );
case 'revweibull'
% REVERSE WEIBULL WITH EXPONENT K
    linkfun = revweibull_link( K, lims );
case 'weibull'
% WEIBULL WITH EXPONENT K
    linkfun = weibull_link( K, lims );
otherwise
% CANONICAL (LOGIT) LINK IF UNKNOW FUNCTION    
    warning('unknown link function: changed to canonical')
    linkfun = logit_link( lims );
end

%%%%% FIT
pfit = glmval(b,x,linkfun );