function [A,tol_bet,stop_val,maxiter,force_rank,max_rank,verbose] =  parse_input_demo(A,varargin)
	p = inputParser;
	p.addRequired('A',@issparse);
	p.addParamValue('tol_bet',0.001,@(x) isa(x,'double'));
	p.addParamValue('stop_val',1e-10,@(x) isa(x,'double'));
	p.addParamValue('maxiter',100,@(x) mod(x,1)==0);
	p.addParamValue('force_rank',0,@(x) mod(x,1)==0);		
	p.addParamValue('max_rank',0,@(x) mod(x,1)==0);		
	p.addParamValue('verbose',true,@(x) isa(x,'logical'));	
	p.parse(A,varargin{:}{:});
	
	A = p.Results.A;
	tol_bet = p.Results.tol_bet;
	stop_val = p.Results.stop_val;
	maxiter = p.Results.maxiter;
	force_rank = p.Results.force_rank;
	max_rank = p.Results.max_rank;	
	verbose = p.Results.verbose;	
	
end