function [r] = colley(T,N)
%% Colley Ranking Model
%
% INPUT: T = n-by-n symmetric matrix
% T(i,j)=-number of games played by i and j
% =0, otherwise
% T(i,i)=0
% N = n-by-3 vector
% N(i,1) = number of wins by team i
% N(i,2) = number of loses by team i
% N(i,3) = total number of games played by team i
%
% OUTPUT: r = vector of ratings

teams=size(T,1);
b=ones(teams,1);
r=zeros(teams,1);
% SET UP COLLEY MATRIX C
C=T+2*eye(teams,teams)+diag(N(:,3));
% SET UP VECTOR b
b=b+0.5*(N(:,1)-N(:,2));
r=inv(C)*b;