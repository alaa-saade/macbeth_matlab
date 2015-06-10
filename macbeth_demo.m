function RMSE = macbeth_demo(varargin)   	
	% MaCBetH : Matrix Completion with the Bethe Hessian
    %
	% Companion paper: Matrix Completion from Fewer Entries: Spectral Detectability 
	% and Rank Estimation.
    % 
	% This is a demo code. run `macbeth_demo;` to run an example with default parameters.
	% macbeth_demo accepts keyword arguments to set a subset of the parameters.
	% For instance, to manually set epsilon to 15 and the rank to 4, 
	% run macbeth_demo('epsilon',15,'rank',4).
	% 
	% List of available keywords 
    %
	% *n* : number of lines of the random matrix to be completed (default 500)
	% *m* : number of columns of the random matrix to be completed (default 500)
	% *rank* : rank of the matrix to be completed (default 2)
	% *epsilon* : average number of revealed entries per row or column (see paper) (default 20)
	% *Delta* : variance of gaussian additive noise (default 0)
	% *stop_val* : stoping value of the minFunc solver. The optimization will stop 
	% if the RMSE on the observed entries is smaller than stop_val (default 1e-10)
	% *maxiter* : (approximate) maximum number of iterations of the minFunc (default 150)
	% *tol_bet* : tolerance of the numerical solver for the parameter beta (default 1e-4)
	% *force_rank* : set to true to force the algorithm to use the correct rank. 
	% By default, force_rank is set to false and Macbeth tries to infer the correct rank. 
	% *max_rank* : maximum possible rank when inferring the rank (default rank+1)
	% *verbose* : set to false to prevent the code from talking (default true)

	addpath('./subroutines');
	[n,m,rank,epsilon,Delta,stop_val,maxiter,tol_bet,force_rank,verbose,max_rank] =  parse_input_demo(varargin);

	if max_rank == 0
		max_rank = rank+1;
	end

	X = randn(n,rank);
	Y = randn(m,rank);

	[obs,noise] = gen_obs(n,m,epsilon,Delta);
	true_A = X*Y';
	A_obs = true_A.*obs + noise;
	
	if verbose
		str = sprintf('Observation matrix computed : rank %d, %dx%d matrix,',rank,n,m);
		disp(str);
		str = sprintf('with %1.2f%% observed entries (%d)', nnz(A_obs)/(n*m)*100,nnz(A_obs));
		disp(str);
	end

	if force_rank == false
		[X_inferred,Y_inferred,r] = complete(A_obs,'tol_bet',tol_bet,'stop_val',stop_val,'maxiter',maxiter,'max_rank',max_rank,'verbose',verbose);
	else 
		[X_inferred,Y_inferred,r] = complete(A_obs,'tol_bet',tol_bet,'stop_val',stop_val,'maxiter',maxiter,'force_rank',rank,'verbose',verbose);
	end
	if r == 0
		RMSE = sqrt(mean(mean((true_A - A_obs).^2)));
	else
		RMSE = sqrt(mean(mean((true_A - X_inferred*Y_inferred').^2)));
	end
	if verbose
		str = sprintf('Reconstruction error (RMSE) on full matrix with a rank %d approach : %1.2e',r,RMSE);
		disp(str);
	end
end

function [obs,noise] = gen_obs(n,m,epsilon,Delta)
	% generates sparsity pattern of observed entries
	% and noise
	density = epsilon/sqrt(n*m);
	idx = randperm(m*n,round(density*m*n));     
	[i,j] = ind2sub([n m],idx);                 
	obs = sparse(i,j,1,n,m);
	if Delta > 1
		noise_vec = sqrt(Delta)*randn(size(i));
		noise = sparse(i,j,noise_vec,n,m);
	else
		noise = sparse(n,m);
	end
end
