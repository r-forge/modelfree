function [b,guess] = binom_g(X,Y,LINK,p,K,initval)
%
% THIS IS AN INTERNAL FUNCTION: USE BINOM_LIMS FOR BEST RESULTS
%
% Maximum likelihood estimates of the parameters of psychometric function
% with guessing rate (GLM). The estimated parameters for the linear part
% are in vector b and the estimated guessing rate is guess
%
% INPUT
%
% X - stimulus levels
% Y - response; two column matrix with the first column corresponding to
% number of successes and the second to number of trials 
%
% OPTIONAL INPUT
% 
% LINK - link function
% p - degree of the polynomial to be fitted on the linear scale; default is
% p=1
% K - Power parameter in Weibull and reverse Weibull models; default is K = 2
% initval - initial value for guess; default is .01
%
% OUTPUT
% 
% b - vector of estimated coefficients for the linear part
% guess - estimated guessing rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROGRAM

%%%% CHECK INPUT PARAMETERS
% First 2 paramaters are mandatory
if (nargin<2)
    error('Design points (X) and results (Y) are mandatory');
end

%%%%
%%%% DEFAULTS
if (nargin<3)
    LINK = 'logit';
    disp('default link function is ''logit''');
end

if (nargin<4)
    p = 1;
    disp('degree of the polynomial to be fitted on the linear scale is 1');
end

if (nargin<5)
    K = 2;
    disp('initial value for K (power parameter in Weibull and reverse Weibull models) is 2');
end

if (nargin<6)
    initval = .01;
    disp('default initial value for guessing rate is 0.01');
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
clear data;
data(1).content = X;
data(2).content = Y;
checkinput( 'psychometricdata', data );
clear data;
checkinput( 'linkfunction', LINK );
checkinput( 'degreepolynomial', p );
% Check that initval is a positive scalar
if length( initval ) > 1 | initval < 0 | initval >= 1
    error( 'guessing rate must be a scalar between 0 and 1' );
end

tmp1(:,1) = X;
X = tmp1;
LX = length(X);
x = repmat(X,1,p).^repmat((1:p),LX,1);

sizeY = size(Y);
if (min(sizeY)==1),
    tmp1(:,1) = Y;
    Y = tmp1;
elseif (min(sizeY)==sizeY(1)),
    Y = Y';
end

if (LX~=max(sizeY)), error('varaible X and Y must have the same length'); end

initval = log(initval/(1-initval));

% GLM ESTIMATION
guess = fminsearch(@(guess) likfun(guess,x,Y,LINK,K),initval,optimset('MaxFunEvals',...
    500,'MaxIter',500,'TolX',1e-3,'TolFun',1e-3));

guess = 1./(1+exp(-guess));

% assign link function
lims = [ guess, 1 ];
switch LINK
case'logit' 
% LOGISTIC
    link = logit_link( lims );
case'probit' 
% PROBIT
    link = probit_link( lims );
case'loglog'
% LOG-LOG
    link = loglog_link( lims );
case'comploglog' 
% COMPLEMENTARY LOG-LOG
    link = comploglog_link( lims );
case 'revweibull'
% REVERSE WEIBULL WITH EXPONENT K
    link = revweibull_link( K, lims );
case 'weibull'
% WEIBULL WITH EXPONENT K
    link = weibull_link( K, lims );
otherwise
% CANONICAL (LOGIT) LINK IF UNKNOW FUNCTION    
    warning('unknown link function: changed to canonical')
    link = logit_link( lims );
end

% GLM
b = glmfit(x,Y,'binomial','link',link);
    
% % % % % % % % % % % % % % % % % % 
% % % INTERNAL FUNCTIONS % % % % % 
% % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%% LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = likfun(guess,x,Y,LINK,K)

guess = 1./(1+exp(-guess));

% assign link function
lims = [ guess, 1 ];
switch LINK
case'logit' 
% LOGISTIC
    link = logit_link( lims );
case'probit' 
% PROBIT
    link = probit_link( lims );
case'loglog'
% LOG-LOG
    link = loglog_link( lims );
case'comploglog' 
% COMPLEMENTARY LOG-LOG
    link = comploglog_link( lims );
case 'revweibull'
% REVERSE WEIBULL WITH EXPONENT K
    link = revweibull_link( K, lims );
case 'weibull'
% WEIBULL WITH EXPONENT K
    link = weibull_link( K, lims );
otherwise
% CANONICAL (LOGIT) LINK IF UNKNOW FUNCTION    
    warning('unknown link function: changed to canonical')
    link = logit_link( lims );
end

% GLM
b = glmfit(x,Y,'binomial','link',link);

% FITTED PROBABILITIES
fitted = glmval(b,x,link);
fitted(fitted<=guess) = guess + eps;
fitted(fitted>=1) = 1 - eps;

% LIKELIHOOD
res = -(Y(:,1)' * log(fitted) + (Y(:,2) - Y(:,1))' * log(1 - fitted));

%%%%%%% END LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%