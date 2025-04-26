# Functions to run PDE simulations of the hydrodynamic equations of the active particle lattice gas

using TensorOperations, LinearAlgebra

## utility 
function midpoint_bond_diff_1d(f::Vector{Float64}; Nx::Int64 = 100, Lx::Float64 = 1.0)

	grad_f::Vector{Float64} = zeros(Nx)

	for x₁ in 1:Nx
		## 1 direction
		y₁::Int64  = (x₁ + Nx)%Nx + 1
		grad_f[x₁] = Nx*(f[y₁] - f[x₁])/Lx
	end
	return grad_f
end

function midpoint_bond_diff_2d(f::Matrix{Float64}; Nx::Int64 = 100, Lx::Float64 = 1.0)

	grad_f::Matrix{Float64} = zeros(Nx, 3)

	for x₁ in 1:Nx
		## 1 direction
		y₁::Int64 = (x₁ + Nx)%Nx + 1
		grad_f[x₁, :] = Nx*(f[y₁, :] - f[x₁, :])/Lx
	end
	return grad_f
end

function midpoint_bond_av_1d(f::Vector{Float64}; Nx::Int64 = 100)
	av_f::Vector{Float64} = zeros(Nx)

	for x₁ in 1:Nx
		## 1 direction
		y₁::Int64 = (x₁ + Nx)%Nx + 1
		av_f[x₁] = (f[y₁] + f[x₁])/2
	end
	return av_f
end

function site_div_1d(f::Vector{Float64}; Nx::Int64 = 100, Lx::Float64 = 1.0)

	div_f::Vector{Float64} = zeros(Nx)

	for x₁ in 1:Nx
		## 1 direction
		y₁::Int64 = (x₁ + Nx)%Nx + 1
		div_f[y₁] += Nx*(f[y₁] - f[x₁])/Lx
	end
	return div_f
end

function site_div_2d(f::Matrix{Float64}; Nx::Int64 = 100, Lx::Float64 = 1.0)

	div_f::Matrix{Float64} = zeros(Nx, 3)

	for x₁ in 1:Nx
		## 1 direction
		y₁::Int64 = (x₁ + Nx)%Nx + 1
		div_f[y₁, :] += Nx*(f[y₁, :] - f[x₁, :])/Lx
	end

	return div_f
end

function flip_term(f::Matrix{Float64}; Nx::Int64 = 100)
	flip_term::Matrix{Float64} = zeros(Nx, 3)
	flip_term[:, 1] = f[:, 1] - f[:, 2]
	flip_term[:, 2] = f[:, 2] - f[:, 1]
	return flip_term
end

function p(x::Float64; logtol = 1e-10, γ = 0.0)
	if x < 0.0
		x = logtol
	elseif x ≥ 1.0
		x = 1.0 - logtol
	end
	c2::Float64 = sqrt((π-2)/((π-1)*(26+(-11+π)*π)))
	c3::Float64 = -(-4+π)sqrt(1+8*(π-3)/(26+(-11+π)*π))/2
	p1::Float64 = (1-π+2*(-3+π)*x)
	p2::Float64 = -2 + 2*π - (-2+π)*(-1+π)*x + (-3+π)*(-2+π)*x^2
	return c3*log(-1-c2*p1) - c3*log(1-c2*p1) + (1 - π)*log(1-x) + 0.5*(-2+π)*log(p2)
end

function mob(f::Matrix{Float64}, ρ::Vector{Float64})
	ds::Vector{Float64} = self_diff.(ρ)
	return f .* ds
end

function upwind(U::Float64, mb_down::Float64, mb_up::Float64)
	return (U > 0.0 ? U .* mb_down : U .* mb_up)
end
#

# coefficients
function self_diff(ρ::Float64; logtol::Float64 = 1e-10)
	α::Float64 = π/2 - 1;
	if ρ ≤ 0.0
		ρ = logtol
	elseif ρ>1.0
		ρ = 1.0
	end
	return (1-ρ) .* (α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ + 1)
end

function self_diff_prime(ρ::Float64; logtol::Float64 = 1e-10)
	α::Float64 = π/2 - 1;
	if ρ ≤ 0.0
		ρ = logtol
	elseif ρ>1.0
		ρ = 1.0
	end
	return - (α*(2*α-1)/(2*α+1)*ρ .^ 2 - α*ρ .+ 1) + (-ρ .+ 1)*(2*α*(2*α-1)/(2*α+1)*ρ - α);
end

ds(x) = self_diff(x)
dsp(x) = self_diff_prime(x)

function mag(f::Matrix{Float64})
	return f[:, 2] - f[:, 1]
end

function coeff_s(rho::Float64, ds::Float64)
	return ((rho*ds)>0 ? (1-rho-ds)/(rho*ds) : 0)
end

function coeff_mag_s(f::Matrix{Float64}, ρ::Vector{Float64})
	m::Vector{Float64}     = mag(f);
	ds::Vector{Float64}    = self_diff.(ρ);
	s::Vector{Float64}     = coeff_s.(ρ, ds);
	mag_s::Vector{Float64} = s .* m
	return mag_s
end
#

## running pde 
function new_pde_param(DT::Float64, v0::Float64, DR::Float64, Δx::Float64, Lx::Float64, ϕa::Float64, ϕp::Float64, δt::Float64, δ::Float64; T::Float64 = 0.001, name::String = "test", save_interval::Float64 = 0.001, save_on::Bool = false)
	param::Dict{String, Any} = Dict{String, Any}()
	Nx::Int64 = Int64(Lx/Δx ÷ 1)
	@pack! param = DT, v0, DR, Δx, Lx, ϕa, ϕp, δt, δ, T, name, Nx, save_interval, save_on
	return param
end

function initiate_uniform_pde(ϕa::Float64, ϕp::Float64, Nx::Int64 = 100)
	f::Matrix{Float64} = zeros(Nx, 3)
	f[:, 1:2] = fill(ϕa/2, (Nx, 2))
	f[:, 3] = fill(ϕp, (Nx))
	return f
end

function U_velocities(f::Matrix{Float64}, ρ::Vector{Float64}; Nx::Int64 = 100, Lx::Float64 = 1.0, DT::Float64 = 1.0, v0::Float64 = 10.0)
	logtol::Float64 = log(1e-10);

	eθ::Array{Float64, 2} = reshape([-1 1 0], 1, 3)

	logmf::Matrix{Float64} = map(x -> (x>0 ? log(x) : logtol), f);
	p_rho::Vector{Float64} = p.(ρ) #functon p is labelled W in the pdf

	U::Matrix{Float64} = -DT*midpoint_bond_diff_2d(logmf .+ p_rho; Nx = Nx, Lx = Lx) .+ v0*midpoint_bond_av_1d(coeff_mag_s(f, ρ); Nx = Nx) .+ v0*eθ

	return U
end

function F_fluxes(U::Matrix{Float64}, moba::Matrix{Float64}; Nx::Int64 = 100)
	F::Matrix{Float64} = zeros(Nx, 3);

	for x₁ in 1:Nx
		local y₁
		## 1 direction
		y₁::Int64 = (x₁ + Nx)%Nx + 1
		F[x₁, :] = upwind.(U[x₁, :], moba[x₁, :], moba[y₁, :])
	end
	return F
end

function time_step!(t::Float64, f::Matrix{Float64}; δt::Float64 = δt, Nx::Int64 = 100, Lx::Float64 = 1.0, DT::Float64 = 1.0, v0::Float64 = 10.0, DR::Float64 = 1.0)
	ρ::Vector{Float64} = sum(f; dims = 2)[:, 1];

	U::Matrix{Float64}    = U_velocities(f, ρ; Nx = Nx, Lx = Lx, DT = DT, v0 = v0);
	mobf::Matrix{Float64} = mob(f, ρ);
	F::Matrix{Float64}    = F_fluxes(U, mobf; Nx = Nx);

	a::Float64 = maximum(abs.(U));

	tempu::Float64 = 1/(6*a*Nx);
	dt::Float64 = min(δt, tempu);

	f -= dt*(site_div_2d(F; Nx = Nx, Lx = Lx) + DR*flip_term(f; Nx = Nx))
	t += dt
	return t, f
end
#

## saving funcitons 
function pde_save_name(param::Dict{String, Any}, t::Float64)
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, name, save_interval = param
	s = round(t; digits = Int64(-log10(save_interval) ÷ 1))
	return datadir("pm_pdes_pro", "$(name)/[DT,v0,DR,Δx,Lx,ϕa,ϕp]=$([DT,v0,DR,Δx,Lx,ϕa,ϕp])/t=$(s).jld2")
end

function speed_save_name(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, name, save_interval = param
	return datadir("pm_pdes_pro", "speeds", "$(name)/[DT,v0,DR,Δx,Lx,ϕa,ϕp]=$([DT,v0,DR,Δx,Lx,ϕa,ϕp]).jld2")
end

function pde_time_series_save_name(param::Dict{String, Any}, t::Float64)
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, name, save_interval = param
	s = round(t; digits = Int64(-log10(save_interval) ÷ 1))
	return datadir("pm_pdes_pro", "$(name)/[DT,v0,DR,Δx,Lx,ϕa,ϕp]=$([DT,v0,DR,Δx,Lx,ϕa,ϕp])/T=$(s)_Δt=$(save_interval).jld2")
end

function load_pde(param::Dict{String, Any}, t::Float64)
	filename::String = pde_save_name(param::Dict{String, Any}, t::Float64)
	data::Dict{String, Any} = load(filename)
	@unpack t, f = data
	return t, f
end

function load_compress_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, save_interval, save_on = param
	t_saves::Vector{Float64} = []
	f_saves::Vector{Matrix{Float64}} = []
	data::Dict{String, Any} = Dict()

	try
		try
			filename::String = pde_time_series_save_name(param, T)
			println(filename)
			data = load(filename)
			@unpack t_saves, f_saves = data
			println(filename)
			println("fast load pde")
			return t_saves, f_saves  # Return immediately after successful load
		catch
			filename::String = pde_time_series_save_name(param, T-save_interval)
			println(filename)
			data = load(filename)
			@unpack t_saves, f_saves = data
			println(filename)
			println("fast load pde (fallback)")
			return t_saves, f_saves  # Return immediately after successful load
		end
	catch
		println("full load pde")
		s = 0.0
		t = 0.0
		while s<T
			try
				t, f = load_pde(param, s)
				push!(f_saves, f)
				push!(t_saves, t)
				s += save_interval
			catch
				s += save_interval
			end
		end
		if t > 0.0
			filename::String = pde_time_series_save_name(param, t)
			data = Dict("f_saves" => f_saves, "t_saves" => t_saves)
			safesave(filename, data)
			println("saved")
		end
	end

	return t_saves, f_saves
end

function silent_load_compress_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, save_interval, save_on = param
	t_saves::Vector{Float64} = []
	f_saves::Vector{Matrix{Float64}} = []
	data::Dict{String, Any} = Dict()

	try
		try
			filename::String = pde_time_series_save_name(param, T)
			data = load(filename)
		catch
			filename::String = pde_time_series_save_name(param, T-save_interval)
			data = load(filename)
		end
		@unpack t_saves, f_saves = data
		# println("fast load")
	catch
		# println("full load")
		s = 0.0
		t = 0.0
		while s<T
			try
				t, f = load_pde(param, s)
				push!(f_saves, f)
				push!(t_saves, t)
				s += save_interval
			catch
				s += save_interval
			end
		end
		if t > 0.0
			filename::String = pde_time_series_save_name(param, t)
			data = Dict("f_saves" => f_saves, "t_saves" => t_saves)
			safesave(filename, data)
			# println("saved")
		end
	end

	return t_saves, f_saves
end

function pde_vid_save_name(param::Dict{String, Any}, t::Float64)
	@unpack DT, v0, DR, N, Lx, Ly, ϕa, ϕp, name, save_interval = param
	s = round(t; digits = Int64(-log10(save_interval) ÷ 1))
	filename = plotsdir("vids", "pm_pdes_vids", "$(name)/[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])/t=$(s).mp4")
	pathname = plotsdir("vids", "pm_pdes_vids", "$(name)/[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])")
	return filename, pathname
end

function load_last_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
	t::Float64 = 0.0;
	s::Float64 = T;
	loaded::Bool = false

	while s>0.0
		try
			t, f = load_pde(param, s)
			loaded = true
			s = -1.0
			println("loaded at t=$(t)")
		catch
			loaded = false
			s += -save_interval
		end
	end

	return loaded, f, t
end

function quiet_load_last_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
	t::Float64 = 0.0;
	s::Float64 = T;
	loaded::Bool = false

	while s>0.0
		try
			t, f = load_pde(param, s)
			loaded = true
			s = -1.0
		catch
			loaded = false
			s += -save_interval
		end
	end

	return loaded, f, t
end
#

## pertubation
function lin_pert_values(param; wave_num = 1, wave_choice = 3)
	@unpack DT, v0, DR, Lx, ϕa, ϕp = param
	ω = 2*π*wave_num/Lx;
	Pe = v0;
	ϕ = ϕa + ϕp;
	ϕ0 = 1 - ϕ;
	ds = self_diff(ϕ);
	dsp = self_diff_prime(ϕ);
	DD = (1-ds)/ϕ
	s = DD - 1
	W = [ -ω^2             0          -im*ω*Pe*ϕ0;
		-ω^2*ϕa*DD      -ω^2*ds     -im*ω*Pe*(ϕa*s+ds);
		-im*ω*Pe*ϕa*dsp -im*ω*Pe*ds -ω^2*ds-2]
	values, vectors = eigen(W)
	return ω, values[wave_choice], vectors[:, wave_choice]
end

function dist_from_unif(f, param)
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt, δ = param
	return sqrt(sum((f[:, 1] .- ϕa/2) .^ 2 + (f[:, 2] .- ϕa/2) .^ 2 + (f[:, 3] .- ϕp) .^ 2)/Nx)
end

function dist_from_unifs(f_saves, param)
	return [dist_from_unif(f, param) for f in f_saves]
end

function node_sym(f)
	return maximum(f - f[end:-1:1, [2, 1, 3]])
end

function perturb_pde!(f::Matrix{Float64}, param::Dict{String, Any}; wave_num = 1)
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt, δ, pert = param

	# Create initial perturbation
	if pert == "rand"
		pertf = 2*rand(Nx, 3) .- 1
	elseif pert == "double"
		ω, value, vector = lin_pert_values(param; wave_num = wave_num, wave_choice = 3)

		wave  = exp.((1:Nx)*(wave_num*im*2*π/Nx))
		pertf = zeros(Nx, 3)

		pertf[:, 1] = real.(wave*(vector[2] - vector[3])/2)
		pertf[:, 2] = real.(wave*(vector[2] + vector[3])/2)
		pertf[:, 3] = real.(wave*(vector[1]-vector[2]))

		pertf += pertf[end:-1:1, [2, 1, 3]]
		print("sym pert: $(node_sym(pertf))")

	else
		ω, value, vector = lin_pert_values(param; wave_num = wave_num)

		wave = exp.((1:Nx)*(wave_num*im*2*π/Nx))
		pertf = zeros(Nx, 3)

		pertf[:, 1] = real.(wave*(vector[2] - vector[3])/2)
		pertf[:, 2] = real.(wave*(vector[2] + vector[3])/2)
		pertf[:, 3] = real.(wave*(vector[1]-vector[2]))

		println("max pert")
	end

	# 1. Normalize by sqrt(Nx) to make scaling independent of grid size
	c = norm(pertf)/sqrt(Nx)
	pertf = δ*pertf/c

	# 2. Ensure zero integral for active and passive components
	# For a function discretized on Nx points over length Lx, 
	# the average value equals the integral divided by Lx
	active_mean = sum(pertf[:, 1] + pertf[:, 2])/Nx
	passive_mean = sum(pertf[:, 3])/Nx

	# Subtract means to ensure zero integral
	pertf[:, 1] .-= active_mean/2  # Distribute correction equally between active components
	pertf[:, 2] .-= active_mean/2
	pertf[:, 3] .-= passive_mean

	# 3. Ensure physical bounds for all components
	# Base densities (assuming equal split of active density)
	ϕa1 = ϕa/2
	ϕa2 = ϕa/2

	# Calculate scaling factors needed for each constraint
	sf_values = []

	# Check individual active component 1 bounds
	max_violation_a1_upper = maximum(ϕa1 .+ pertf[:, 1]) - 1.0
	max_violation_a1_lower = -minimum(ϕa1 .+ pertf[:, 1])

	if max_violation_a1_upper > 0
		push!(sf_values, (1.0 - ϕa1) / maximum(pertf[:, 1]))
	end

	if max_violation_a1_lower > 0
		push!(sf_values, ϕa1 / maximum(-pertf[:, 1]))
	end


	# Check individual active component 2 bounds
	max_violation_a2_upper = maximum(ϕa2 .+ pertf[:, 2]) - 1.0
	max_violation_a2_lower = -minimum(ϕa2 .+ pertf[:, 2])

	if max_violation_a2_upper > 0
		push!(sf_values, (1.0 - ϕa2) / maximum(pertf[:, 2]))
	end

	if max_violation_a2_lower > 0
		push!(sf_values, ϕa2 / maximum(-pertf[:, 2]))
	end


	# Check total active component bounds
	max_violation_a_upper = maximum(ϕa .+ pertf[:, 1] .+ pertf[:, 2]) - 1.0
	max_violation_a_lower = -minimum(ϕa .+ pertf[:, 1] .+ pertf[:, 2])

	if max_violation_a_upper > 0
		push!(sf_values, (1.0 - ϕa) / maximum(pertf[:, 1] + pertf[:, 2]))
	end
	if max_violation_a_lower > 0
		push!(sf_values, ϕa / maximum(-(pertf[:, 1] + pertf[:, 2])))
	end

	# Check passive component bounds
	max_violation_p_upper = maximum(ϕp .+ pertf[:, 3]) - 1.0
	max_violation_p_lower = -minimum(ϕp .+ pertf[:, 3])

	if max_violation_p_upper > 0
		push!(sf_values, (1.0 - ϕp) / maximum(pertf[:, 3]))
	end
	if max_violation_p_lower > 0
		push!(sf_values, ϕp / maximum(-pertf[:, 3]))
	end

	# Check total density bounds
	total_pert = sum(pertf, dims = 2)
	max_violation_total_upper = maximum(ϕa .+ ϕp .+ total_pert) - 1.0
	max_violation_total_lower = -minimum(ϕa .+ ϕp .+ total_pert)

	if max_violation_total_upper > 0
		push!(sf_values, (1.0 - (ϕa + ϕp)) / maximum(total_pert))
	end
	if max_violation_total_lower > 0
		push!(sf_values, (ϕa + ϕp) / maximum(-total_pert))
	end

	# Apply the most restrictive scaling factor if any constraint is violated
	if !isempty(sf_values)
		sf = minimum(sf_values) * (1.0 - 1e-5)  # Add a small buffer for numerical stability
		pertf *= sf
	end

	# Verify zero integral constraint is still satisfied after scaling
	# (scaling shouldn't affect this, but let's be extra careful)
	active_mean = sum(pertf[:, 1] + pertf[:, 2])/Nx
	passive_mean = sum(pertf[:, 3])/Nx

	pertf[:, 1] .-= active_mean/2
	pertf[:, 2] .-= active_mean/2
	pertf[:, 3] .-= passive_mean

	# Add perturbation to original field
	f += pertf

	return f
end

## solving functions

function run_new_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
	f = perturb_pde!(f, param);

	t::Float64 = 0.0;
	s::Float64 = save_interval

	#inital save
	if save_on
		filename::String        = pde_save_name(param, t)
		data::Dict{String, Any} = Dict("f" => f, "t" => t)
		safesave(filename, data)
	end

	while t < T
		while t < s
			t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
		end
		#save snapshot
		if save_on
			filename = pde_save_name(param, t)
			data     = Dict("f" => f, "t" => t)
			safesave(filename, data)
		end
		s += save_interval
	end
	return t, f
end

function run_new_pde_sym(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
	f = perturb_pde!(f, param);
	t::Float64 = 0.0;
	s::Float64 = save_interval
	f = (f + f[end:-1:1, :]) ./ 2

	#inital save
	if save_on
		filename::String        = pde_save_name(param, t)
		data::Dict{String, Any} = Dict("f" => f, "t" => t)
		safesave(filename, data)
	end

	while t < T
		while t < s
			t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
			f = (f + f[end:-1:1, :]) ./ 2
		end
		#save snapshot
		if save_on
			filename = pde_save_name(param, t)
			data     = Dict("f" => f, "t" => t)
			safesave(filename, data)
		end
		s += save_interval
	end
	return t, f
end

function run_current_pde(param::Dict{String, Any}, dt::Float64, f::Matrix{Float64}, t::Float64)
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	s::Float64 = t + save_interval
	t_end::Float64 = t+dt

	#inital save
	if save_on
		filename::String = pde_save_name(param, t)
		data::Dict{String, Any} = Dict("f" => f, "t" => t)
		safesave(filename, data)
	end

	while t < t_end
		while t < s
			t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
		end
		#save snapshot
		if save_on
			filename = pde_save_name(param, t)
			data     = Dict("f" => f, "t" => t)
			safesave(filename, data)
		end
		s += save_interval
	end
	return t, f
end

function run_current_pde_sym(param::Dict{String, Any}, dt::Float64, f::Matrix{Float64}, t::Float64)
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	s::Float64 = t + save_interval
	t_end::Float64 = t+dt

	#inital save
	if save_on
		filename::String = pde_save_name(param, t)
		data::Dict{String, Any} = Dict("f" => f, "t" => t)
		safesave(filename, data)
	end

	while t < t_end
		while t < s
			t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
			f = (f + f[end:-1:1, :]) ./ 2
		end
		#save snapshot
		if save_on
			filename = pde_save_name(param, t)
			data     = Dict("f" => f, "t" => t)
			safesave(filename, data)
		end
		s += save_interval
	end
	return t, f
end

function load_and_run_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
	t::Float64 = 0.0;
	s::Float64 = T;
	loaded::Bool = false

	while s>=0.0
		try
			t, f = load_pde(param, s)
			loaded = true
			s = -1.0
		catch
			loaded = false
			s += -save_interval
		end
	end

	if loaded
		println("load at t = $(t)")
		s = t + save_interval
		while t < T
			while t < s
				t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
			end
			#save snapshot
			if save_on
				filename = pde_save_name(param, t)
				data     = Dict("f" => f, "t" => t)
				safesave(filename, data)
			end
			s += save_interval
		end
	else
		println("all loading failed; running new pde")
		t, f = run_new_pde(param)
	end
	return t, f
end

function run_and_pert_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
	f = perturb_pde!(f, param);
	t::Float64 = 0.0;
	s::Float64 = save_interval

	#inital save
	if save_on
		filename::String        = pde_save_name(param, t)
		data::Dict{String, Any} = Dict("f" => f, "t" => t)
		safesave(filename, data)
	end

	pert_count = 1
	δ = 0.01
	pert = "rand"
	@pack! param = δ, pert
	while t < T
		while t < s
			t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
		end
		if pert_count > 100
			f = perturb_pde!(f, param)
			pert_count = 1
		else
			pert_count += 1
		end
		#save snapshot
		if save_on
			filename = pde_save_name(param, t)
			data     = Dict("f" => f, "t" => t)
			safesave(filename, data)
		end
		s += save_interval
	end
	return t, f
end

function force_load_and_pert_pde(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	# configuration
	f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
	t::Float64 = 0.0;
	s::Float64 = save_interval;
	loaded::Bool = false

	try
		t, f = load_pde(param, 0.0)
		loaded = true
	catch
		loaded = false
	end

	pert_count = 1
	δ = 0.01
	pert = "rand"
	@pack! param = δ, pert
	if loaded
		while t < T
			while t < s
				t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
			end
			if pert_count > 100
				f = perturb_pde!(f, param)
				pert_count = 1
			else
				pert_count += 1
			end
			#save snapshot
			if save_on
				filename = pde_save_name(param, t)
				data     = Dict("f" => f, "t" => t)
				safesave(filename, data)
			end
			s += save_interval
		end
	else
		println("all loading failed; abort")
	end
	return t, f
end
#

## stretch fns 
function double_sol(param, f)
	@unpack Nx, Lx, Δx = param
	NNx = Int64(2*Nx)
	g = zeros(NNx, 3)
	g[2:2:NNx, :] = f
	g[1:2:NNx, :] = (f + circshift(f, (1, 0))) / 2
	Δx = Δx/2
	Nx = NNx
	@pack! param = Δx, Nx
	return g, param
end

function relax_sol(param, f, t; threshold = 1e-4)
	@unpack Lx, save_interval = param
	dc = threshold/Lx + 1
	while abs(dc) > threshold/Lx
		t, f = run_current_pde(param, save_interval, f, t)
		_, _, dc = f_dot(param, f)
	end
	return f, t
end

function get_stretch_param(Lx)
	param = get_grid_param(21, 11)
	@unpack Nx = param
	param["save_interval"] = 100.0
	param["name"] = "soliton_stretch"
	param["Lx"] = Float64(Lx)
	param["Δx"] = Float64(Lx/Nx)
	return param
end

function densify(Lx, ΔX; save_interval = 1.0, threshold = 1e-6)

	param = get_stretch_param(Lx)
	@pack! param = save_interval

	loaded, f, t = quiet_load_last_pde(param)
	t += 2000.0

	while param["Δx"] > ΔX
		f, param = double_sol(param, f)
		println("relaxing Δx = $(param["Δx"])")
		f, t = relax_sol(param, f, t; threshold = threshold)
	end

	filename = pde_save_name(param, t)
	data     = Dict("f" => f, "t" => t)
	safesave(filename, data)

	return f, t
end

function get_dense_param(Lx, ΔX; save_interval = 1.0)
	param = get_stretch_param(Lx)
	@pack! param = save_interval

	while param["Δx"] > ΔX
		param = double_param(param)
	end
	return param
end

function double_param(param)
	@unpack Nx, Lx, Δx = param
	NNx = Int64(2*Nx)
	Δx = Δx/2
	Nx = NNx
	@pack! param = Δx, Nx
	return param
end
#
