function [f,F,g,G,Gp,loglik]=LDSsmooth(v,A,B,CovH,CovV,CovP,meanP,meanH,meanV,varargin)
%LDSSMOOTH Linear Dynamical System : Filtering and Smoothing
%[f,F,g,G,Gp]=LDSsmooth(v,A,B,CovH,CovV,CovP,meanP,meanH,meanV,<cellout>)
%
% Inputs:
% v : V x T matrix of observations
% A : transition matrix
% B : emission matrix
% CovH : transition covariance
% CovV : emission covariance
% CovP : initial (prior) covariance
% meanP : initial (prior) mean
% meanH : mean of the hidden transition
% meanV : mean of the observation
%
% Outputs:
% f,F : filtered mean and covariance
% g,G : smoothed mean and covariance
% Gp : smoothed cross moment <h_t h_{t+1}|v_{1:T}> t=1..T-1
% loglik : log likelihood of the sequence v
% By default the outputs are stored in arrays, eg F(:,:,t). Calling
% with cellout set to 1 returns cells eg F{t}
cellout=0; if nargin==10; cellout=varargin{1}; end
[f2,F2,loglik]=LDSforward(v,A,B,CovH,CovV,meanH,meanV,CovP,meanP);
[g2,G2,Gp2]=LDSbackward(v,A,B,f2,F2,CovH,meanH);
[V T]=size(v); H=size(CovH,1);
if cellout
    f=mat2cell(f2,H,ones(1,T));
    F=mat2cell(F2,H,H,ones(1,T));
    g=mat2cell(g2,H,ones(1,T));
    G=mat2cell(G2,H,H,ones(1,T));
    Gp=mat2cell(Gp2,H,H,ones(1,T-1));    
else
    f=f2; F=F2; g=g2; G=G2; Gp=Gp2;
end