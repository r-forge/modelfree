function b = binomfit_lims(r, m, x, p, link, lims, K)
%
% Fits a binomial generalised linear model with lower and upper constraints
% as in lims
%
%INPUT
%
% r - number of successes in points x
% m - number of trials in points x 
% x - stimulus levels
%
% OPTIONAL INPUT
%
% p - degree of the polynomial to be fitted on the linear scale; default is
% p=1 
% link - link function; choose from: 'logit', 'probit', 'loglog',
% 'comploglog', 'weibull', 'revweibull'; default is the canonical link 'logit' 
% lims - two column vector specifying guessing and 1-lapsing rates; default
% is [0,1] 
% K - the exponent for weibull and revweibull link function; default is 2
%
% OUTPUT
%
% b - estimated coefficients for the linear part; column vector of length
% p+1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PROGRAM

%%%% CHECK INPUT PARAMETERS
% First 2 paramaters are mandatory
if (nargin<3)
    error('Check input. First 3 arguments are mandatory');
end

%%%%
%%%% DEFAULTS
if (nargin<4)
    p = 1;
    disp('degree of the polynomial to be fitted on the linear scale is 1');
end

if (nargin<5)
    link = 'logit';
    disp('default link function is ''logit''');
end

if (nargin<6)
    lims = [0 1];
    disp('default lower and upper limits are 0 and 1');
end

if (nargin<7)
    K = 2;
    if strcmp(link, 'weibull')
        disp('default exponent for Weibull link function is 2');
    elseif strcmp(link, 'revweibull')
        disp('default exponent for reverse Weibull link function is 2');
    end
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
clear data;
data(1).content = x;
data(2).content = r;
data(3).content = m;
checkinput( 'psychometricdata2', data );
clear data;
checkinput( 'degreepolynomial', p );
checkinput( 'linkfunction', link );
checkinput( 'guessingandlapsing', lims );
checkinput( 'exponentk', K );

%%%%
%%%% INITIALS VALUES

% column vectors
tmp1(:,1) = x;
x = tmp1;

Lx = length(x);

% create matrix with all powers (1 to p) of X (for use in GLMFIT)
X = repmat(x,1,p).^repmat((1:p),Lx,1);

%%%%%%
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

%  fit the GLM model
b = glmfit( X, [r m], 'binomial', 'link', linkfun);