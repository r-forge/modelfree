function [ b, lims ] = binom_lims( r, m, x, gl, link, p, K, initval )
%
% Maximum likelihood estimates of the parameters of the psychometric function
% with guessing and lapsing rates or only guessing rate. The
% estimated parameters for the linear part are in vector b and the
% estimated limits are in lims.
%
% INPUT
%
% r - number of successes in points x
% m - number of trials in points x 
% x - stimulus levels
%
% OPTIONAL INPUT
% 
% gl - Indicator, calulate only guessing if 'guess' and both guessing and
% lapsing otherwise
% link - link function
% p - degree of the polynomial to be fitted on the linear scale; default is
% p=1 
% K - Power parameter in Weibull and reverse Weibull models; default is K=2
% initval - initial value for limits; default is [.01 .99] if guessing and
% lapsing rates are calculated and .01 if only guessing rate is calculated
%
% OUTPUT
% 
% b - vector of estimated coefficients for the linear part
% lims - estimated guessing and 1-lapsing rates, or only the guessing rate (if
% gl = 'guess')
%
% Created by Ivan Marin-Franch, 20/11/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROGRAM

%%%% CHECK INPUT PARAMETERS

% First 3 paramaters are mandatory
if (nargin<3)
    error('Check input. First 3 arguments are mandatory');
end

%%%%
%%%% DEFAULTS
if ( nargin<4 )
    gl = 'both';
    disp('both guessing and lapsing rates are calculated');
end

if ~strcmp( gl, 'both' )  && ...
   ~strcmp( gl, 'guess' )
    error( 'wrong value for guessing/lapsing indicator gl' );
end

if (nargin<5)
    link = 'logit';
    disp('default link function is ''logit''');
end

if (nargin<6)
    p = 1;
    disp('degree of the polynomial to be fitted on the linear scale is 1');
end

if (nargin<7)
    K = 2;
    disp('initial value for K (power parameter in Weibull and reverse Weibull models) is 2');
end

if (nargin<8)
    if strcmp( gl, 'both' )
        initval = [.01 .99];
        disp('default initial values for guessing and 1 - lapsing rates are 0.01 and 0.99');
    else
        initval = .01;
        disp('default initial value for guessing rate is 0.01');
    end
end


%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
clear data;
data(1).content = x;
data(2).content = r;
data(3).content = m;
checkinput( 'psychometricdata2', data );
clear data;
checkinput( 'linkfunction', link );
checkinput( 'degreepolynomial', p );

% Check initial values for guessing or guessing/lapsing and call internal
% functions binom_gl or binom_g depending on indicator gl
if strcmp( gl, 'both' )
    checkinput( 'guessingandlapsing', initval );
    [ b, lims ] = binom_gl( x, [r m], link, p, K, initval );
else
% Check that initval is a positive scalar
    if length( initval ) > 1 || initval <= 0 || initval >= 1
        error( 'guessing rate must be a scalar between 0 and 1' );
    end
    [ b, lims ] = binom_g( x, [r m], link, p, K, initval );
end