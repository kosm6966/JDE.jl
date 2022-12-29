"""

`run_optim(P0,LB,UB; NT,NS,RT,RUNTIME,VERB)` runs the Simulated Annealing
program to optimize the parameters given in P0 by recursively minimizing χ².
Inputs P0, LB, UB are required. Results are printed to screen.

P0: Array{Float64,1}. Initial Values for Parameters. [N,khop,KHT]
LB: Array{Float64,1}. Lower Bounds on Parameters. [lb1,lb2,lb3]
UB: Array{Float64,1}. Upper Bounds on Parameters. [ub1,ub2,ub3]
NT: Int64. Reduce temperature every NT*NS*dim(P0) evaluations.
NS: Int64. Adjust bounds every NS*dim(P0) evaluations.
RT: Float64. Cooling Rate. When temp changes, new temp is T_new = RT * T_old
RUNTIME: Int64. Terminate run after RUNTIME seconds.

"""
function run_optim(P0, LB, UB; NT=10, NS=20, RT=0.85, RUNTIME=30)
	test = TestData()
	b0 = test.x
	measured = test.prompt
	Ninv = 1. / length(measured)

	experiment = Experiment(distribution="./test/testdata/powder_distribution.csv",field=test.x)
	spec = Spectrum(basis=:Adiabatic)
	S = Vector{Float64}(undef,length(b0))

	function objective_function(P0)
		dimer = BridgedDimer(d=P0[1],x=P0[2],β=P0[3],g=P0[4])
		S .= quintetspectrum(dimer,spec,experiment)

		residual =  sum( (measured .- S*P0[5]).^2 ) * Ninv
		return residual
	end

	optimize(objective_function, LB, UB, P0,
			SAMIN(;
				nt = NT,     # 10; reduce temperature every NT*NS*dim(P0) evaluations
				ns = NS,     # 20; adjust bounds every NS*dim(P0) evaluations
				rt = RT,     # 0.85; geometric temperature reduction factor: when temp changes, new temp is T=RT*T0
				neps = 5,    # number of previous best values the final result is compared to
				f_tol = 1e-3, # 1e-3, the required tolerance level for function value comparisons
				x_tol = 1e-2, # 1e-2, the required tolerance level for x
				coverage_ok = false, # if false, increase temperature until initial parameter space is covered
				verbosity = 3), # scalar: 0, 1, 2 or 3 (default = 0). If 0 - don't print to screen. If 3 - print everything to screen.
			Optim.Options(
				iterations = 10^6,
				  time_limit = RUNTIME, # seconds; 6 days: 518400
				show_trace = true,
				store_trace = false,
				extended_trace = true))

	println("\n run complete \n")

	return nothing
end

function makerandp0(LB,UB,n)
	out = zeros(n,length(LB))
	for i in eachindex(LB)
		out[:,i] = LB[i] .+ (UB[i] - LB[i]) .*rand(n)
	end
        out .= round.(out,digits=2)
        return out
end