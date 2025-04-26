# particle simulations of the active particle lattice gas
using StatsBase, Random

## utility 
function weight_index(i::Int64, j::Int64, k::Int64; Nx::Int64 = 1, Ny::Int64 = 1)
	return i + Nx*(j-1) + Nx*Ny*(k-1)
end
#

## running simulation 
function new_sim_param(DT::Float64, v0::Float64, DR::Float64, N::Int64, Lx::Float64, Ly::Float64, ϕa::Float64, ϕp::Float64; T::Float64 = 0.001, name::String = "test", save_interval::Float64 = 0.001, save_on::Bool = false)
	param::Dict{String, Any} = Dict{String, Any}()
	N₁::Int64, N₂::Int64 = Int64(Lx*N ÷ 1), Int64(Ly*N ÷ 1)
	@pack! param = DT, v0, DR, N, Lx, Ly, ϕa, ϕp, T, name, N₁, N₂, save_interval, save_on
	return param
end

function get_jump_rate(η::Array{Float64, 3}, x₁::Int64, x₂::Int64, y₁::Int64, y₂::Int64, jump::Int64, DT::Float64, v0::Float64, N::Int64, N₁::Int64, N₂::Int64)
	is_valid::Int64 = η[x₁, x₂, 1]*(1-η[y₁, y₂, 1])
	bias::Float64 = 0.0
	if jump == 2 #left
		bias = -η[x₁, x₂, 2]
	elseif jump == 3 #right
		bias = η[x₁, x₂, 2]
	end
	return is_valid*(DT*N^2 + v0*bias*N/2)
end

function initiate_uniform(ϕa::Float64, ϕp::Float64, N₁, N₂)
	η::Array{Float64, 3} = zeros(N₁, N₂, 2)
	# x₁, x₂, spin = i, j , k
	w::Weights{Float64, Float64, Vector{Float64}} = weights([ϕa/2, ϕa/2, ϕp, (1 - ϕa - ϕp)])
	particle::Vector{Vector{Int64}} = [[1, 1], [1, -1], [1, 0], [0, 0]]
	for x₁ in 1:N₁, x₂ in 1:N₂
		η[x₁, x₂, :] = sample(particle, w)
	end
	return η
end

function initiate_weights(η::Array{Float64, 3}, N₁::Int64, N₂::Int64, DT::Float64, v0::Float64, DR::Float64, N::Int64)
	w::Array{Float64, 3} = zeros(N₁, N₂, 5)
	jumps::Vector{Tuple{Int64, Int64}} = [(i, j) for i in -1:1:1, j in -1:1:1 if i^2+j^2 == 1];
	for x₁ in 1:N₁, x₂ in 1:N₂, jump in 1:4
		local y₁::Int64, y₂::Int64
		# find adjacent site
		y₁, y₂ = ((x₁, x₂) .+ jumps[jump] .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
		w[x₁, x₂, jump] = get_jump_rate(η, x₁, x₂, y₁, y₂, jump, DT, v0, N, N₁, N₂)
	end
	for x₁ in 1:N₁, x₂ in 1:N₂
		local jump::Int64
		jump = 5
		w[x₁, x₂, jump] = DR*η[x₁, x₂, 2]^2
	end
	return weights(w)
end

function model_step!(η::Array{Float64, 3}, w::Weights{Float64, Float64, Vector{Float64}}, t::Float64, N₁::Int64, N₂::Int64, DT::Float64, v0::Float64, DR::Float64, N::Int64, jumps::Vector{Tuple{Int64, Int64}})
	#update total propensity
	prop::Float64 = sum(w)
	#update time
	t += randexp()/prop
	#select jump
	x::Array{Tuple{Int64, Int64, Int64}, 3} = [(i, j, k) for i in 1:N₁, j in 1:N₂, k in 1:5]
	x₁::Int64, x₂::Int64, jump::Int64 = sample(x, w)
	if jump == 5
		#execute jump
		η[x₁, x₂, 2] = - η[x₁, x₂, 2]
		#update rates
		for jump in 1:4
			local y₁::Int64, y₂::Int64
			# find adjacent site
			y₁, y₂ = ((x₁, x₂) .+ jumps[jump] .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
			w[weight_index(x₁, x₂, jump; Nx = N₁, Ny = N₂)] = get_jump_rate(η, x₁, x₂, y₁, y₂, jump, DT, v0, N, N₁, N₂)
		end
	else
		local y₁::Int64, y₂::Int64
		# find adjacent site
		y₁, y₂ = ((x₁, x₂) .+ jumps[jump] .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
		# swap particles
		η[x₁, x₂, :], η[y₁, y₂, :] = η[y₁, y₂, :], η[x₁, x₂, :]
		# set hop rates
		for jump in 1:4
			local z₁::Int64, z₂::Int64
			#update adjacent sites to x
			z₁, z₂ = ((x₁, x₂) .+ jumps[jump] .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
			w[weight_index(x₁, x₂, jump; Nx = N₁, Ny = N₂)] = 0.0 # as x is empty
			w[weight_index(z₁, z₂, 5-jump; Nx = N₁, Ny = N₂)] = get_jump_rate(η, z₁, z₂, x₁, x₂, 5-jump, DT, v0, N, N₁, N₂)

			#update adjacent sites to y
			z₁, z₂ = ((y₁, y₂) .+ jumps[jump] .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
			w[weight_index(y₁, y₂, jump; Nx = N₁, Ny = N₂)] = get_jump_rate(η, y₁, y₂, z₁, z₂, jump, DT, v0, N, N₁, N₂)
			w[weight_index(z₁, z₂, 5-jump; Nx = N₁, Ny = N₂)] = get_jump_rate(η, z₁, z₂, y₁, y₂, 5-jump, DT, v0, N, N₁, N₂)
		end
		# set flip rates
		w[weight_index(x₁, x₂, 5; Nx = N₁, Ny = N₂)] = 0.0
		w[weight_index(y₁, y₂, 5; Nx = N₁, Ny = N₂)] = DR*η[y₁, y₂, 2]^2
	end
	return t
end
#

## saving funcitons 
function sim_save_name(param::Dict{String, Any}, t::Float64)
	@unpack DT, v0, DR, N, Lx, Ly, ϕa, ϕp, name, save_interval = param
	s = round(t; digits = Int64(-log10(save_interval) ÷ 1))
	return datadir("pm_sims_raw", "$(name)", "[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])", "t=$(s).jld2")
end

function sim_time_series_save_name(param::Dict{String, Any}, t::Float64)
	@unpack DT, v0, DR, N, Lx, Ly, ϕa, ϕp, name, save_interval = param
	s = round(t; digits = Int64(-log10(save_interval) ÷ 1))
	return datadir("pm_sims_pro", "$(name)", "[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])", "T=$(s)_Δt=$(save_interval).jld2")
end

function load_sim(param::Dict{String, Any}, t::Float64)
	filename::String = sim_save_name(param::Dict{String, Any}, t::Float64)
	data::Dict{String, Any} = load(filename)
	@unpack η, t = data
	return t, η
end

function load_compress_sim(param::Dict{String, Any})
	@unpack DT, v0, DR, N, Lx, ϕa, ϕp, T, name, N₁, N₂, save_interval, save_on = param
	t_saves::Vector{Float64} = []
	η_saves::Vector{Array{Float64, 3}} = []
	data::Dict{String, Any} = Dict()

	try
		try
			filename::String = sim_time_series_save_name(param, T)
			println(filename)
			data = load(filename)
		catch
			filename::String = sim_time_series_save_name(param, T-save_interval)
			println(filename)
			data = load(filename)
		end
		@unpack t_saves, η_saves = data
		println("fast load sim")
	catch
		println("full load sim")
		s = 0.0
		t = 0.0
		while s<T
			try
				t, η = load_sim(param, s)
				push!(η_saves, η)
				push!(t_saves, t)
				s += save_interval
			catch
				s += save_interval
			end
		end
		if t > 0.0
			filename::String = sim_time_series_save_name(param, t)
			data = Dict("η_saves" => η_saves, "t_saves" => t_saves)
			safesave(filename, data)
			println("saved")
		end
	end

	return t_saves, η_saves
end
#

## simulation functions 

function run_new_sim(param::Dict{String, Any})
	@unpack DT, v0, DR, N, Lx, ϕa, ϕp, T, name, N₁, N₂, save_interval, save_on = param
	# configuration
	η::Array{Float64, 3} = initiate_uniform(ϕa, ϕp, N₁, N₂);
	w::Weights{Float64, Float64, Vector{Float64}} = initiate_weights(η, N₁, N₂, DT, v0, DR, N);
	t::Float64 = 0.0;
	s::Float64 = save_interval
	jumps::Vector{Tuple{Int64, Int64}} = [(i, j) for i in -1:1:1, j in -1:1:1 if i^2+j^2 == 1];

	#inital save
	if save_on
		filename::String        = sim_save_name(param, t)
		data::Dict{String, Any} = Dict("η" => η, "t" => t)
		safesave(filename, data)
	end

	while t < T
		while t < s
			t = model_step!(η, w, t, N₁, N₂, DT, v0, DR, N, jumps);
		end
		#save snapshot
		if save_on
			filename = sim_save_name(param, t)
			data     = Dict("η" => η, "t" => t)
			safesave(filename, data)
		end
		s += save_interval
	end
	return t, η, w
end

function run_current_sim(param::Dict{String, Any}, dt::Float64, η::Array{Float64, 3}, w::Weights{Float64, Float64, Vector{Float64}}, t::Float64)
	@unpack DT, v0, DR, N, Lx, ϕa, ϕp, name, N₁, N₂, save_interval, save_on = param
	# configuration
	s::Float64 = t + save_interval
	jumps::Vector{Tuple{Int64, Int64}} = [(i, j) for i in -1:1:1, j in -1:1:1 if i^2+j^2 == 1];
	t_end::Float64 = t+dt

	#inital save
	if save_on
		filename::String = sim_save_name(param, t)
		data::Dict{String, Any} = Dict("η" => η, "t" => t)
		safesave(filename, data)
	end

	while t < t_end
		while t < s
			t = model_step!(η, w, t, N₁, N₂, DT, v0, DR, N, jumps);
		end
		#save snapshot
		if save_on
			filename = sim_save_name(param, t)
			data = Dict("η" => η, "t" => t)
			safesave(filename, data)
		end
		s += save_interval
	end
	return t, η, w
end

function load_and_run_sim(param::Dict{String, Any})
	@unpack DT, v0, DR, N, Lx, ϕa, ϕp, T, name, N₁, N₂, save_interval, save_on = param
	# configuration
	η::Array{Float64, 3} = initiate_uniform(ϕa, ϕp, N₁, N₂);
	w::Weights{Float64, Float64, Vector{Float64}} = initiate_weights(η, N₁, N₂, DT, v0, DR, N);
	t::Float64 = 0.0;
	s::Float64 = T;
	loaded::Bool = false
	jumps::Vector{Tuple{Int64, Int64}} = [(i, j) for i in -1:1:1, j in -1:1:1 if i^2+j^2 == 1];

	while s>0.0
		try
			t, η = load_sim(param, s)
			w = initiate_weights(η, N₁, N₂, DT, v0, DR, N);
			loaded = true
			s = -1.0
		catch
			loaded = false
			#println("load failed at t = $(s)")
			s += -save_interval
		end
	end

	if loaded
		println("load at t = $(t)")
		s = t + save_interval
		while t < T
			while t < s
				t = model_step!(η, w, t, N₁, N₂, DT, v0, DR, N, jumps);
			end
			#save snapshot
			if save_on
				filename = sim_save_name(param, t)
				data     = Dict("η" => η, "t" => t)
				safesave(filename, data)
			end
			s += save_interval
		end
	else
		println("all loading failed; running new simulation")
		t, η, w = run_new_sim(param)
	end
	return t, η, w
end
#

## proccessing funcitons 

function local_average(η, ϵ, N, N₁, N₂)
	f::Array{Float64, 3} = zeros(N₁, N₂, 3)
	jumps::Vector{Tuple{Int64, Int64}} = [(i, j) for i in (-N₁):1:N₁, j in (-N₂):1:N₂ if i^2+j^2 ≤ ϵ*N];
	njumps::Int64 = length(jumps)
	for x₁ in 1:N₁, x₂ in 1:N₂
		for jump in jumps
			local y₁::Int64, y₂::Int64
			# find adjacent site
			y₁, y₂ = ((x₁, x₂) .+ jump .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
			f[x₁, x₂, 1] += (η[y₁, y₂, 2]^2-η[y₁, y₂, 2])/2
			f[x₁, x₂, 2] += (η[y₁, y₂, 2]^2+η[y₁, y₂, 2])/2
			f[x₁, x₂, 3] += (η[y₁, y₂, 1]^2-η[y₁, y₂, 2]^2)
		end
	end
	return f/njumps
end

function local_average_1d(η, ϵ, N, N₁, N₂)
	f::Matirx{Float64} = zeros(N₁, 3)
	jumps::Vector{Tuple{Int64, Int64}} = [(i, j) for i in (-N₁):1:N₁, j in (-N₂):1:N₂ if i^2+j^2 ≤ ϵ*N];
	njumps::Int64 = length(jumps)
	for x₁ in 1:N₁, x₂ in 1:N₂
		for jump in jumps
			local y₁::Int64, y₂::Int64
			# find adjacent site
			y₁, y₂ = ((x₁, x₂) .+ jump .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
			f[x₁, 1] += (η[y₁, y₂, 2]^2-η[y₁, y₂, 2])/2
			f[x₁, 2] += (η[y₁, y₂, 2]^2+η[y₁, y₂, 2])/2
			f[x₁, 3] += (η[y₁, y₂, 1]^2-η[y₁, y₂, 2]^2)
		end
	end
	return f/(njumps*N₂)
end

function local_average_timeseries(η_saves, ϵ, N, N₁, N₂)
	num_t = length(η_saves)
	ft::Array{Float64, 3} = zeros(num_t, N₁, 3)
	jumps::Vector{Tuple{Int64, Int64}} = [(i, j) for i in (-N₁):1:N₁, j in (-N₂):1:N₂ if i^2+j^2 ≤ ϵ*N];
	njumps::Int64 = length(jumps)
	for (i, η) in enumerate(η_saves)
		for x₁ in 1:N₁, x₂ in 1:N₂
			for jump in jumps
				local y₁::Int64, y₂::Int64
				# find adjacent site
				y₁, y₂ = ((x₁, x₂) .+ jump .+ (N₁-1, N₂-1)) .% (N₁, N₂) .+ (1, 1)
				ft[i, x₁, 1] += (η[y₁, y₂, 2]^2-η[y₁, y₂, 2])/2
				ft[i, x₁, 2] += (η[y₁, y₂, 2]^2+η[y₁, y₂, 2])/2
				ft[i, x₁, 3] += η[y₁, y₂, 1]^2-η[y₁, y₂, 2]^2
			end
		end
	end
	return ft ./ (njumps*N₂)
end

#
