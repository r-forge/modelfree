function checkinput(type, data)
%
% THIS IS AN INTERNAL FUNCTION ONLY
%
% Pool of routines that check robustness of input of parameters for all
% functions of modelfree package.
%
%INPUT
%
% type - Type of checking
% data - Input data to be checked
%
% Created by Ivan Marin-Franch, 13/11/2008

switch type
    case 'designpoints'
        checkdesignpoints( data );
    case 'psychometricdata'
        checkpsychometricdata( data(1).content, data(2).content );
    case 'psychometricdata2'
        checkpsychometricdata2( data(1).content, data(2).content, data(3).content );
    case 'degreepolynomial'
        checkdegreepolynomial( data );
    case 'linkfunction'
        checklinkfunction( data );
    case 'guessingandlapsing'
        checkguessingandlapsing( data );
    case 'bootstrapreplications'
        checkbootstrapreplications( data );
    case 'exponentk'
        checkexponentk( data );
    case 'minmaxbandwidth'
        checkminmaxbandwidth( data );
    case 'bandwidth'
        checkbandwidth( data );
    case 'kernel'
        checkkernel( data );
    case 'maxiter'
        checkmaxiter( data );
    case 'tolerance'
        checktolerance( data );
    case 'method'
        checkmethod( data );
end

end

function checkdesignpoints( x )

[ rx, cx ] = size( x );
if cx ~= 1
    error('vector of stimulus levels should be column vector');
end

end

function checkpsychometricdata(x,y)

% Check format of x and y
[ rx, cx ] = size( x );
[ ry, cy ] = size( y );

if rx < 2 || ry < 2
    error( 'Minimum number of points is 2' );
end

if rx ~= ry
    error( 'Dimension mismatch. Number of design points must equal number of reponses' );
end

if cx ~= 1
    error( 'Vector of design points should be column vector' );
end

if cy ~= 2
    error( 'Response matrix must have two columns: number of successes and number of trials' );
end

if any( y(:,1) > y(:,2) )
    error( 'Wrong response matrix. Number of successes cannot be larger than number of trials' );
end

end

function checkpsychometricdata2( xdes, k, m )

% Check format of x and y
[ rx, cx ] = size( xdes );
[ rk, ck ] = size( k );
[ rm, cm ] = size( m );

if rx ~= rk || rx ~= rm || rk ~= rm
    error( 'The number of design points for design points, successes and trials must be the same' );
end

if rx < 2 || rk < 2 || rm < 2
    error( 'Minimum number of points is 2' );
end

if ( rx ~= rk ) || ( rx ~= rm )
    error( 'Dimension mismatch. Number of design points must equal number of reponses' );
end

if cx ~= 1
    error( 'Vector of design points should be column vector' );
end

if ck ~= 1
    error( 'Vector of number of successes should be column vector' );
end

if cm ~= 1
    error( 'Vector of number of trials should be column vector' );
end

if any( k > m )
    error( 'Wrong response matrix. Number of successes cannot be larger than number of trials' );
end

end

function checkdegreepolynomial( p )

if p <= 0 | round( p ) ~= p | length( p ) > 1
    error( 'degree of polynomial has to be a positive integer' );
end

end

function checklinkfunction( LINK )

if ~strcmp( LINK, 'logit' )      && ...
   ~strcmp( LINK, 'probit' )     && ...
   ~strcmp( LINK, 'loglog' )     && ...
   ~strcmp( LINK, 'comploglog' ) && ...
   ~strcmp( LINK, 'weibull' )    && ...
   ~strcmp( LINK, 'revweibull' )
    error( '%s is not an allowed link function', char( LINK ) );
end

end

function checkguessingandlapsing( lims )

[ rlims clims ] = size( lims );

if rlims ~= 1 || clims ~= 2
    error( 'lims has be a row vector with two values: Lower and Upper limits' );
end
if any( lims < 0 ) || any( lims > 1 )
    error( 'guessing and lapsing rates must be between 0 and 1' );
end
if lims(1) >= lims(2)
    error( 'guessing rate cannot be larger than nor equal to lapsing' );
end

end

function checkbootstrapreplications( N )

if N <= 0 | round( N ) ~= N | length( N ) > 1
    error( 'number of bootstrap replications has to be a positive integer' );
end

end

function checkexponentk( k )

if ( length( k ) > 1 || k <= 0 )
    error( 'exponent for Weibull or reverse Weibull link function must be a positive scalar' );
end

end

function checkminmaxbandwidth( H )

[ rH cH ] = size( H );

if rH ~= 1 || cH ~= 2
    error( 'H has be a row vector with two values: Lower and Upper limits' );
end
if H(1) >= H(2)
    error( 'minimum bandwidth cannot be larger than nor equal to the upper one' );
end

if H(1) <= 0
    error( 'bandwidth cannot be negative nor zero' );
end

end

function checkbandwidth( h )

if length( h ) > 1 | h <= 0
    error( 'bandwidth must be a positive scalar' );
end

end

function checkkernel( ker )

if ~ischar( ker )
    ker = func2str( ker );
end

if ~strcmp( ker, 'normpdf' )      && ...
   ~strcmp( ker, 'epanechnikov' ) && ...
   ~strcmp( ker, 'triangular' )   && ...
   ~strcmp( ker, 'tricube' )      && ...
   ~strcmp( ker, 'bisquare' )     && ...
   ~strcmp( ker, 'uniform' )
    error( '"%s" is not an allowed kernel', ker );
end

end
    
function checkmaxiter( maxiter )

if maxiter <= 0 | round( maxiter ) ~= maxiter | length( maxiter ) > 1
    error( 'maximum number of iterations has to be a positive integer' );
end

end

function checktolerance( tol )

if ( length( tol ) > 1 || tol <= 0 )
    error( 'tolerance level must be a positive scalar' );
end

end

function checkmethod( method )

[ rm cm ] = size( method );

if rm > 1
    error( 'Choose only one method or ''all'' to calculate bandwidth' );
end

if ~strcmp( method, 'ISE' )        && ...
   ~strcmp( method, 'ISEeta' )     && ...
   ~strcmp( method, 'deviance' ) && ...
   ~strcmp( method, 'all' )
    error( '"%s" is a wrong loss function', method );
end

end