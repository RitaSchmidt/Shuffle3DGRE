function [Mask] = EllipticMask(Length1, Length2)

  % get k-space "indices" centered around zero.
  I0 = Ordering.Utils.FFTCenterIndex(Length1) ;
  J0 = Ordering.Utils.FFTCenterIndex(Length2) ;
  [ii, jj]=ndgrid((1:Length1) - I0, ...
                  (1:Length2) - J0) ;

  % Scaling so we can handle rectangular matrices. Similar to 'a', 'b' in
  % defintions of ellipses, only their reciprocals (e.g., a --> 1/a), I
  % think. 
  a = Length1/min(Length1, Length2) ;
  b = Length2/min(Length1, Length2) ;

  % Effective radius of all points ("elliptical"), acounting for the
  % rectangular shape of the matrix.
  rEff = sqrt(b^2*ii.^2 + a^2*jj.^2) ;

  % Actually Define the mask
  % We would like to use floor below (for odd lengths), but slight
  % rounding errors might throw us off. Since the results should be
  % either integer or half integer, subtractinbg 0.1 and rounding, should
  % work as well. (Assuming rounding errors are much smaller than 0.1)
  % Mask = (rEff <= floor(max(Length1/2, Length2/2))) ;
  RLim = round(max(Length1/2, Length2/2) -0.1) ;
  Mask = (rEff <= RLim + eps(RLim)) ;


end