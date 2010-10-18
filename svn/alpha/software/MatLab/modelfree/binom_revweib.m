function [b,K] = binom_revweib( r, m, x, p, initK, lims)
%
% Maximum likelihood estimates of the parameters of the reverse Weibull model
% for the psychometric function. The estimated parameters for the linear
% part are in vector b and the estimated exponent is K.
%
%
% INPUT
%
% r - number of successes in points x
% m - number of trials in points x 
% x - stimulus levels
%
% OPTIONAL INPUT
% 
% p - degree of the polynomial to be fitted on the linear scale; default is
% p=1 
% initK - initial value for K (power parameter in reverse Weibull model);
% default is initK=2 
% lims - two column vector with guessing and 1-lapsing rate; default is
% [0,1] 
%
% OUTPUT
% 
% b - vector of estimated coefficients for the linear part
% K - estimate of the power parameter in the reverse Weibull model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % PROGRAM

%%%% CHECK INPUT PARAMETERS
% First 2 paramaters are mandatory
if (nargin<3)
    error('Check input. First 3 arguments are mandatory');
end

%%%% DEFAULTS
if (nargin<4)
    p = 1;
    disp('degree of the polynomial to be fitted on the linear scale is 1');
end

if (nargin<5)
    initK = 2;
    disp('initial value for K (power parameter in reverse Weibull model) is 2');
end

if (nargin<6)
    lims = [0 1];
    disp('default lower and upper limits are 0 and 1');
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
clear data;
data(1).content = x;
data(2).content = r;
data(3).content = m;
checkinput( 'psychometricdata2', data );
clear data;
checkinput( 'degreepolynomial', p );
checkinput( 'guessingandlapsing', lims );
checkinput( 'exponentk', initK );

tmp1(:,1) = x;
x = tmp1;
Lx = length(x);
X = repmat(x,1,p).^repmat((1:p),Lx,1);

% GLM ESTIMATION

K = fminsearch(@(K) likfun(K,X,[r m],lims),initK,optimset('MaxFunEvals',5000,...
    'MaxIter',5000,'TolX',1e-3,'TolFun',1e-3));
K = .05+exp(K);

%%%%%%%%%
%%% SET LINK
link = revweibull_link( K, lims );

%%%%%%
% GLM
b = glmfit( X, [r m], 'binomial', 'link', link);

% % % % % % % % % % % % % % % % % % 
% % % INTERNAL FUNCTIONS % % % % % 
% % % % % % % % % % % % % % % % % % 


%%%%%%%%%%%%%%%%%% LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = likfun(K,x,Y,lims)

K = .05+exp(K);

%%% SET LINK
link = revweibull_link( K, lims );

% GLM
b = glmfit(x,Y,'binomial','link',link);

% FITTED PROBABILITIES
fitted = glmval(b,x,link);
fitted(fitted<=lims(1)) = lims(1) + eps;
fitted(fitted>=lims(2)) = lims(2) - eps;

% LIKELIHOOD
res = -(Y(:,1)' * log(fitted) + (Y(:,2) - Y(:,1))' * log(1 - fitted));

%%%%%%% END LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%