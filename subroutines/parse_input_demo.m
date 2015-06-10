function [n,m,rank,epsilon,Delta,stop_val,maxiter,tol_bet,force_rank,verbose,max_rank] =  parse_input_demo(varargin)
	p = inputParser;
	p.addParameter('n',500,@(x) mod(x,1)==0);
	p.addParameter('m',500,@(x) mod(x,1)==0);
	p.addParameter('rank',2,@(x) mod(x,1)==0);
	p.addParameter('epsilon',20,@(x) isa(x,'double'));
	p.addParameter('Delta',0,@(x) isa(x,'double'));
	p.addParameter('stop_val',1e-10,@(x) isa(x,'double'));
	p.addParameter('maxiter',150,@(x) mod(x,1)==0);
	p.addParameter('tol_bet',0.0001,@(x) isa(x,'double'));
	p.addParameter('force_rank',false,@(x) isa(x,'logical'));		
	p.addParameter('verbose',true,@(x) isa(x,'logical'));		
	p.addParameter('max_rank',0,@(x) mod(x,1)==0);
	p.parse(varargin{:}{:});
	
	n = p.Results.n;
	m = p.Results.m;
	rank = p.Results.rank;
	epsilon = p.Results.epsilon;
	Delta = p.Results.Delta;
	stop_val = p.Results.stop_val;
	maxiter = p.Results.maxiter;
	tol_bet = p.Results.tol_bet;
	force_rank = p.Results.force_rank;
	verbose = p.Results.verbose;
	max_rank = p.Results.max_rank;
end
