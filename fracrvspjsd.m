
%
%
% compute the fractional reverse permutation JS distance between two time series 
%
% 
%
%
% inputs
%
%   x - an input 1D signal, a sequence
%
%   m - the embedding dimension, a positive integer
%   
%   tau - the time delay, a positive integer
%
%   alpha -   -1 <= frac <= 1   fractional exponent
%
%
%
%
% outputs
%
%   frpjsdsc - the score of fractional reverse JS distance
%     
%
%
%
% e.g.
%
%  [ aa ] = fracrvspjsd( rrits, 3, 1, 0.8 );
%
%
%
%
% refs
%
%   [1] Financial time series analysis based on fractional and multiscale permutation entropy
%
%
%
%
%
%
%
% author - koke yao 1064686304@qq.com
% 
%
%
%





function [ frpjsdsc ] = fracrvspjsd( x, m, tau, alpha )



%%  vectrorization


x = x( : ); 




%%  list all possible permutation patterns according to ref [4]


% M = perms( 1 : 1 : m ); % only for check
% % % M = flipud( M );




%%  phase space reconstruction XX


Nx = length( x ); % signal length


if m == 1
    xhkmt = x( 1 : tau : end, 1 );
else
    hkmt = hankel( 1 : tau : Nx-(m-1)*tau, Nx-(m-1)*tau : tau : Nx );
    xhkmt = x( hkmt ); % reconstructed vectors
end

clear  hkmt  x  Nx  


[ ~, xphkmt ] = sort( xhkmt, 2 ); % sort each reconstructed vectors
clear  xhkmt




%%  frequency computation  XX


[ ~, ia, ic ] = unique( xphkmt, 'rows' ); % i - the row indices of the first appearance
xfreq = histcounts( ic, length( ia ) )';
xfreq( length( xfreq ) + 1  : 1 : factorial( m ) ) = 0;
xfreq = xfreq( : );
clear  xphkmt  ia  ic




%%  probability


xp = xfreq ./ sum( xfreq );

yp = ( 1 / factorial( m ) ) * ones( factorial( m ), 1 );

mp = ( xp + yp ) ./ 2;




%%  final defination


frpjsdsc = fracshanen( mp, alpha ) - 0.5 * fracshanen( xp, alpha ) - 0.5 * fracshanen( yp, alpha );




%%  fractinal Shannon entropy


    function  fse = fracshanen( pp, a )


        if  a == 0


            logp = log10( pp );
            logp( isinf( logp ) ) = 0;

            fse = - sum( logp .* pp );


        else


            logp = log10( pp );
            logp( isinf( logp ) ) = 0;

            pna = pp .^ ( - a );
            pna( isinf( pna ) ) = 0;

            fse = - sum( ( ( logp + psi( 1 ) - psi( 1 - a ) ) .* pna ./ gamma( a + 1 ) ) .* pp );


        end

    end



end







