function [ x,y ] = find_node( k1,b1,k2,b2 )
%FIND_NODE Summary of this function goes here
%   Detailed explanation goes here
    if(k1 == k2) 
        error('两条直线平行，没有交点');
    else
        x=(b2-b1)/(k1-k2);
        y=k1*x+b1;
    end
    x=round(x);
    y=round(y);
end

