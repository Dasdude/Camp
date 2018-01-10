function [ x_dbm ] = linear2dbm( x )
%LINEAR2DBM Summary of this function goes here
%   Detailed explanation goes here
x_dbm = 10*log10((x.^2)*1000);

end

