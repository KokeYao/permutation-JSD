
%
%
% compute the permutation JS distance between two time series 
% 
%
%
% inputs
%
%   x - an input 1D signal, a sequence
%
%   y - an input 1D signal, a sequence
%
%   m - the embedding dimension, a positive integer
%   
%   tau - the time delay, a positive integer
%
%
%
% outputs
%
%   npjsdsc - the score of distance
%
%   npjsdsc - the score of normalized distance
%     
%
%
% e.g.
%
%  [ aa ] = pjsd( rrits, qtits, 5, 6 );
%
%
%
% refs
%
%   [1] 
%
%
%
%
% author - koke yao 1064686304@qq.com
% 
%
%
%




function [ pjsdsc, npjsdsc ] = pjsd( x, y, m, tau )



%%  vectrorization

x = x( : ); 
y = y( : ); 




%%  list all possible patterns according to ref [4]

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


xfreq = zeros( factorial( m ), 1 ); % predefine a frequency distribution


for i = 1 : 1 : size( xphkmt, 1 )
    

    mk = xphkmt( i, 1 );
    starting = factorial( m - 1 ) * ( m - mk );
    
    col = 2;
    
    
    while  col < m + 1
                
        mk = xphkmt( i, col );
        vec = sort( xphkmt( i, col : end ), 'descend' ); %
        [ ~ , n ] = find( vec == mk );
        
        starting = starting + factorial( m - col ) * ( n - 1 );
        
        col = col + 1;
        
    end
        
    row = starting + 1;
    xfreq( row, 1 ) = xfreq( row, 1 ) + 1; % add 1 to the frequency distribution
    

    clear  row  starting  col  vec  mk  n
end

clear  i  xphkmt


xfreq = flipud( xfreq );




%%  phase space reconstruction YY


Ny = length( y ); % signal length


if m == 1
    yhkmt = y( 1 : tau : end, 1 );
else
    hkmt = hankel( 1 : tau : Ny-(m-1)*tau, Ny-(m-1)*tau : tau : Ny );
    yhkmt = y( hkmt ); % reconstructed vectors
end

clear  hkmt  y  Ny


[ ~, yphkmt ] = sort( yhkmt, 2 ); % sort each reconstructed vectors
clear  yhkmt




%%  frequency computation  YY


yfreq = zeros( factorial( m ), 1 ); % predefine a frequency distribution


for i = 1 : 1 : size( yphkmt, 1 )
    
    mk = yphkmt( i, 1 );
    starting = factorial( m - 1 ) * ( m - mk );
    
    col = 2;
    
    while  col < m + 1
                
        mk = yphkmt( i, col );
        vec = sort( yphkmt( i, col : end ), 'descend' ); %
        [ ~ , n ] = find( vec == mk );
        
        starting = starting + factorial( m - col ) * ( n - 1 );
        
        col = col + 1;
        
    end
        
    row = starting + 1;
    yfreq( row, 1 ) = yfreq( row, 1 ) + 1; % add 1 to the frequency distribution
    
    clear  row  starting  col  vec  mk
end
clear  i  yphkmt

yfreq = flipud( yfreq );




%%  probability


xp = xfreq ./ sum( xfreq );
yp = yfreq ./ sum( yfreq );
mp = ( xp + yp ) / 2;

xp( ( xp == 0 ) ) = [ ];
yp( ( yp == 0 ) ) = [ ];
mp( ( mp == 0 ) ) = [ ];




%%  final defination


pjsdsc = - sum( mp .* log10( mp ) ) + 0.5 * sum( xp .* log10( xp ) ) + 0.5 * sum( yp .* log10( yp ) );
npjsdsc = pjsdsc / log10( 2 );



end







