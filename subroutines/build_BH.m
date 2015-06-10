function BH = build_BH(A,n,m,tol_bet,c_1,c_2,verbose)
	% lin,col,val = findnz(J)
	% N = length(J(1,:));
	bet = solve_beta(A,tol_bet,c_1,c_2);

	J = [sparse(n,n) A ; A' sparse(m,m)];
	N = length(J(1,:));

	if max(max(abs(tanh(bet*J)))) < 0.99
		if verbose
			disp('Using canonical formulation of the Bethe Hessian');
		end

		BH = -0.5*sinh(2.0*bet*J);
		D = sum(sinh(bet*J).^2);
		D = sparse(1:N,1:N,D,N,N,N);
		BH = BH + speye(N,N) + D;
	else 
		if verbose
			disp('Using signed formulation of the Bethe Hessian for better numerical stability');
		end
		J1 = sign(J);
		c = sqrt(c_1*c_2);
		D = sum(abs(J1),1);
		BH = (c-1)*speye(N)-sqrt(c)*J1+sparse(1:N,1:N,D,N,N);
	end
end

function bet = solve_beta(A,tol,c_1,c_2)
	% lin,col,val = findnz(J)
	bet_min = 0.0;
	bet_max = 1000.0;
	current_value = F( 1/2 * (bet_max + bet_min),A,c_1,c_2);
	
	while bet_max - bet_min > tol

		if current_value > 0
			bet_max = 1/2 * (bet_max + bet_min);
		else 
			bet_min = 1/2 * (bet_max + bet_min);
		end
		current_value = F( 1/2 * (bet_max + bet_min),A,c_1,c_2);
		
	end

	bet = 1/2 * (bet_max + bet_min);
end

function y = F(bet,A,c_1,c_2)
	S = sum(sum(tanh(bet*A).^2));
	num_obs = nnz(A);
	y = sqrt(c_1*c_2)*S/num_obs - 1;
end

