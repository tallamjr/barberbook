function score=BayesNetScore(data,A,nstates,U)
%BAYESNETSCORE Bayes Net score for a dataset data under a Dirichlet prior
% score=BayesNetScore(data,A,nstates,U)
% data is a data matrix with a datum per column
% A is a belief net adjacency matrix
% nstates are the number of states of the variables 1,2,3,...
% U is the hyperprior value for the Dirichlet distribution
% see demoBDscore.m
score=0;
for v=1:size(A,1)
    Pa=parents(A,v);
    score = score+BDscore(data(v,:),data(Pa,:),nstates(v),nstates(Pa),U);
end