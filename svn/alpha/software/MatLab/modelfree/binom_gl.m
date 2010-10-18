function [b,lims] = binom_gl(X,Y,LINK,p,K,initval)
%
% THIS IS AN INTERNAL FUNCTION: USE BINOM_LIMS FOR BEST RESULTS
%
% Maximum likelihood estimates of the parameters of psychometric function
% with guessing and lapsing rates (GLM). The estimated parameters for the
% linear part are in vector b and the estimated limits are in lims.   
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
% p = 1
% K - Power parameter in Weibull and reverse Weibull models; default is K = 2
% initval - initial value for limits; default is [.01 .99] 
%
% OUTPUT
% 
% b - vector of estimated coefficients for the linear part
% lims - estimated guessing and lapsing rates

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
    initval = [.01 .99];
    disp('default initial values for guessing and 1 - lapsing rates are 0.01 and 0.99');
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
clear data;
data(1).content = X;
data(2).content = Y;
checkinput( 'psychometricdata', data );
clear data;
checkinput( 'linkfunction', LINK );
checkinput( 'degreepolynomial', p );
checkinput( 'guessingandlapsing', initval );

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

% GLM ESTIMATION
lims = fminsearch(@(lims) likfun(lims,x,Y,LINK,K),initval,optimset('MaxFunEvals',...
    5000,'MaxIter',5000,'TolX',1e-3,'TolFun',1e-3));

lims = 1./(1+exp(-lims));

% assign link function
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

function res = likfun(lims,x,Y,LINK,K)

lims = 1./(1+exp(-lims));

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
fitted(fitted<=lims(1)) = lims(1) + eps;
fitted(fitted>=lims(2)) = lims(2) - eps;

% LIKELIHOOD
res = -(Y(:,1)' * log(fitted) + (Y(:,2) - Y(:,1))' * log(1 - fitted));

%%%%%%% END LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%