function [CenterIdx] = FFTCenterIndex(N)

% [CenterIdx] = FFTCenterIndex(N)
%
% Returns the "center" index for FFT of length N:
% N odd => CenterIdx = (N+1)/2
% N even => CenterIdx = N/2 + 1
%
% * N positive and non-integer, rounds N first and returns result for
%   roudned value. (Negative case, see below)
% * N = 0 returns 0
% * N < 0 returns minus the result of the matching positive value
% * N can be an array


  Reminder = round(mod(abs(N), 2)) ;
  
  CenterIdx = zeros(size(N)) ;
  
  OddLogicals = (Reminder >= 0.5  & Reminder < 1.5) ;
  
  CenterIdx(OddLogicals) = (abs(N(OddLogicals)) + 1) / 2 ;
  CenterIdx(~OddLogicals) = abs(N(~OddLogicals)) / 2 + 1 ;
  CenterIdx(N==0) = 0 ; % special case of zero
  CenterIdx(N<0) = -CenterIdx(N<0) ; % special case of negative numbers
                                     % Should not happen, but just to be
                                     % safe.
  

  CenterIdx = round(CenterIdx) ; % To be on the safe side

return ;




