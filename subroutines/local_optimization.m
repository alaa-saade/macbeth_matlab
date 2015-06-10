function [inferred_X,inferred_Y] = local_optimization(starting_vec,r,A,n,m,stop_val,maxiter,verbose)
	[I,J,VAL] = find(A);
	Nnz = length(I);

	global count;
	count = 0;

	if (exist('minFunc')==2)    
		disp('Cleaning using minFunc from http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html');
		options.progTol=Nnz*stop_val^2;
		options.maxiter=maxiter;
		options.Display = 'off';
		options.GradObj = 'on';
		options.optTol=-Inf;
		myfun = @(x)costfunc(x,I,J,VAL,n,m,r,verbose);
		minx = minFunc(myfun,starting_vec,options);
	else
		disp('Cleaning using lsqnonlin from the optimization toolbox');
		disp('For much better performance, consider using minFunc from http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html');
		[u,v] = compute_sparsity_pattern_jac(I,J,Nnz,n,m,r);
		myfun = @(x)costfunc_lsqnonlin(x,I,J,VAL,n,m,r,u,v,verbose);
		options = optimoptions('lsqnonlin','Diagnostics','off','Jacobian','on','Display','iter-detailed','DerivativeCheck','off','Algorithm','levenberg-marquardt','ScaleProblem','Jacobian','MaxIter',15,'TolFun',1e-30,'TolX',1e-30); 
		options.Display = 'off';
		minx = lsqnonlin(myfun,starting_vec,[],[],options);
	end
	inferred_X = reshape(minx(1:n*r),n,r);
	inferred_Y = reshape(minx(n*r+1:(n+m)*r),m,r);
end

function [F,Jac] = costfunc_lsqnonlin(x,I,J,VAL,n,m,r,u,v,verbose)
	global count;
	count = count+1;
	Nnz = length(I);
	if nargout == 1
		F = compute_current_estimate_vec(x,I,J,VAL,n,m,r,Nnz);
	else
		[F,JVAL] = compute_current_estimate_vec_wjac(x,I,J,VAL,n,m,r,Nnz);
		Jac = sparse(u,v,JVAL,Nnz,(m+n)*r);
	end

	if verbose 
		rmse = sqrt(1/Nnz*sum(F.^2));
		str = sprintf('Function call #%d,\t RMSE on observed entries = %e',count,rmse);
		disp(str);
	end

end

function [F,grad] = costfunc(x,I,J,VAL,n,m,r,verbose)
	global count;
	count = count+1;
	Nnz = length(I);
	if nargout == 1
		F = compute_current_estimate(x,I,J,VAL,n,m,r,Nnz);
	else
		[F,grad] = compute_current_estimate_wgrad(x,I,J,VAL,n,m,r,Nnz);
	end

	if verbose && mod(count,20) ==1
		str = sprintf('Function call #%d,\t RMSE on observed entries = %e',count,sqrt(1/Nnz*F));
		disp(str);
	end


end


