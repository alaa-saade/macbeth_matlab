function [X0,Y0,r] = infer_BH(BH,n,m,force_rank,max_rank,verbose)

	N = length(BH(1,:));
	Tr = trace(BH);
	BH2 = 1/N*(Tr*speye(N,N) - BH);

	if max_rank >0

		[v,d] = eigs( BH2 , max_rank);
		d = diag(d);
		nconv = length(d);
		[d,id] = sort(d,'descend');
		% d = d[id]
		v = v( :, id );
		r = 0;

		for i=1:max_rank
			if d(i)>1/N*Tr
				r = r + 1;
			end
		end

		r = min(r,nconv);

		if r == 0
			warning('The Bethe Hessian found that the rank is 0 : failure. Will output 0.');
			X0 = 0;
			Y0 = 0;
			return;
		elseif r == max_rank
			w=sprintf('The rank is possibly larger than %d. You might want to increase max_rank',max_rank);
			warning(w);
		end
		if verbose
			str = sprintf('Using an inferred rank equal to %d for the factorization', r);
			disp(str);
		end
		X0 = v(1:n,1:r);
		Y0 = v(n+1:m+n,1:r);
	else
		[v,d] = eigs( BH2,force_rank);
		d = diag(d);
		nconv = length(d);
		[d,id] = sort(d,'descend');
		% d = d[id]
		v = v(:, id);
		r = force_rank;
		if verbose
			str = sprintf('Using the specified rank %d for the factorization',r);
			disp(str);
		end
		X0 = v(1:n,1:r);
		Y0 = v(n+1:m+n,1:r);

	end
end
