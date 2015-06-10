function [inferred_X,inferred_Y,r] = complete(A,varargin)
	% MaCBetH : Matrix Completion with the Bethe Hessian
    %
	% Main completion function. Usage : [inferred_X,inferred_Y,inferred_r] = complete(M)
	% where M is the matrix to be completed. Returns the inferred factors X and Y such 
	% that M is approximately equal to XY', and the inferred rank r. `complete` accepts 
	% keyword arguments to set a subset of the parameters.
	% See macbeth_demo.m for an example on a synthetic low-rank matrix.
	% Note that the input observation matrix should be `centered`.
    %
	% List of available keywords 
    %
	% *tol_bet* : tolerance of the numerical solver for the parameter beta (default 1e-4)
	% *stop_val* : stoping value of the minFunc solver. The optimization will stop 
	% *maxiter* : (approximate) maximum number of iterations of the minFunc (default 150)
	% *force_rank* : set to nonzero value to force Macbeth to use the specified rank. 
	% Either force_rank or max_rank should be set to a nonzero value.
	% *max_rank* : Number of eigenvalues of the hessian to be computed. If all the eigenvalues 
	% computed are negative (i.e. if the inferred rank is larger than max_rank), Macbeth will 
	% give you a warning. Either force_rank or max_rank should be set to a nonzero value.
	% *verbose* : set to false to prevent the code from talking (default true)

[A,tol_bet,stop_val,maxiter,force_rank,max_rank,verbose] = parse_input_complete(A,varargin);

if max_rank==0 && force_rank==0
	error('Either max_rank or force_rank should be >0')
end
if max_rank~=0 && force_rank~=0
	error('Either max_rank or force_rank should be equal to 0')
end
if max_rank<0 || force_rank<0
	error('max_rank and force_rank should be non-negative integer')
end

if max_rank >0 && verbose
	str = strcat('Completion with unspecified rank (max_rank = ',num2str(max_rank),')');
	disp(str);
elseif force_rank>0 && verbose 
	str = strcat('Completion with specified rank (force_rank = ',num2str(force_rank),')');
	disp(str);
end


[n,m]= size(A);
A1 = spones(A);
c_1 = mean(sum(A1,1));
c_2 = mean(sum(A1,2));


BH = build_BH(A,n,m,tol_bet,c_1,c_2,verbose);

if verbose
	disp('Bethe Hessian built');
end

[X0,Y0,r] = infer_BH(BH,n,m,force_rank,max_rank,verbose);
if r == 0
	inferred_X = 0;
	inferred_Y = 0;
else
	if verbose
		disp('Initial inference done, proceeding to local optimization')
	end
	starting_vec = [reshape(X0,n*r,1) ; reshape(Y0,m*r,1)];
	[inferred_X,inferred_Y] = local_optimization(starting_vec,r,A,n,m,stop_val,maxiter,verbose);
end
end