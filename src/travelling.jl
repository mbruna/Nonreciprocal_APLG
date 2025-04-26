# Functions to solve for the travelling solutions of the APLG using the NonlinearSolve package, including the full solution and the inner-outer solutions.

using KernelDensity, KernelDensitySJ, Peaks, Statistics, ForwardDiff;
using NonlinearSolve, DifferentialEquations
using LaTeXStrings

d2(x) = round(x; digits = 2);
d3(x) = round(x; digits = 3);
d4(x) = round(x; digits = 4);
d6(x) = round(x; digits = 6);

## running simulation
function new_param(
	DT::Float64,
	v0::Float64,
	DR::Float64,
	N::Int64,
	Δx::Float64,
	Lx::Float64,
	Ly::Float64,
	ϕa::Float64,
	ϕp::Float64,
	δt::Float64,
	δ::Float64;
	T::Float64 = 0.001,
	name::String = "test",
	pert::String = "lin",
	save_interval::Float64 = 0.001,
	save_on::Bool = false,
)
	param::Dict{String, Any} = Dict{String, Any}()
	N₁::Int64, N₂::Int64 = Int64(Lx*N ÷ 1), Int64(Ly*N ÷ 1)
	Nx::Int64 = Int64(Lx/Δx ÷ 1)
	@pack! param = DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ, T, name, Nx, N₁, N₂, save_interval, save_on, pert
	return param
end

function get_soliton_param(i)
	if i < 12
		DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
		Lx = 20.0
		ϕa = collect(0.37:0.01:0.47)[i]
		ϕp = 0.3
		v0 = 7.5
		T, save_interval, param_name, pert = (Lx^2*5, 10.0, "soliton_speed_check", "lin")
		return new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	elseif i-11 < 11
		DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
		Lx = collect(21.0:1:30.0)[i-11]
		ϕa = 0.37
		ϕp = 0.3
		v0 = 7.5
		Δx = Float64(Lx/400.0)
		T, save_interval, param_name, pert = (Lx^2*5, 10.0, "soliton_speed_check", "lin")
		return new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	elseif i-21 < 16
		DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
		Lx = 20.0
		ϕa = 0.46
		ϕp = collect(0.29:(-0.01):0.15)[i-21]
		v0 = 7.5
		T, save_interval, param_name, pert = (Lx^2*5, 10.0, "soliton_speed_check", "lin")
		return new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	elseif i-36 < 26 # total is 61 for now keep 45 
		DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
		Lx = 20.0
		ϕa = 0.46
		ϕp = 0.3
		v0 = collect(7.6:(0.1):10.0)[i-36]
		T, save_interval, param_name, pert = (Lx^2*5, 10.0, "soliton_speed_check", "lin")
		return new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	else
		print("out of range")
		DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
		Lx = 20.0
		ϕa = 0.46
		ϕp = 0.3
		v0 = 7.5
		T, save_interval, param_name, pert = (Lx^2*5, 10.0, "soliton_speed_check", "lin")
		return new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	end
end

function get_active_param(i)
	DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
	Lx = 20.0
	ϕa = collect(0.3:0.01:0.6)[i]
	ϕp = 0.3
	v0 = 7.5
	T, save_interval, param_name, pert = (2000.0, 10.0, "soliton_vertical_sweep", "lin")
	param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	param["Nx"] = 400
	return param
end

function get_grid_param(i, j)
	DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
	Lx = 20.0
	ϕa = collect(0.3:0.01:0.6)[i]
	ϕp = collect(0.2:0.01:0.4)[j]
	v0 = 7.5
	T, save_interval, param_name, pert = (2000.0, 20.0, "soliton_grid", "lin")
	param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	param["Nx"] = 400
	return param
end

function get_grid_param_wide(i, j)
	DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
	Lx = 20.0
	ϕa = collect(0.2:0.02:0.8)[i]
	ϕp = collect(0.02:0.02:0.4)[j]
	v0 = 7.5
	T, save_interval, param_name, pert = (1000.0, 20.0, "soliton_grid", "lin")
	param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	param["Nx"] = 400
	param["save_interval"] = 200.0
	param["T"] = 1001.0
	param["name"] = "large_pert_grid"
	param["δ"] = 0.1
	param["pert"] = "double"
	return param
end

function get_grid_param_wide(i, j)
	DT, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.01);
	Lx = 20.0
	ϕa = collect(0.2:0.02:0.8)[i]
	ϕp = collect(0.02:0.02:0.4)[j]
	v0 = 7.5
	T, save_interval, param_name, pert = (1000.0, 20.0, "soliton_grid", "lin")
	param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	param["Nx"] = 400
	param["save_interval"] = 200.0
	param["T"] = 1001.0
	param["name"] = "large_pert_grid"
	param["δ"] = 0.1
	param["pert"] = "double"
	return param
end

function cp_param(i)
	DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 25.0, 0.5, 0.1, 0.1, 1e-5, 0.1);
	ϕa, ϕp = [(x, y) for x in 0.1:0.02:1.0, y in 0.0:0.01:0.4 if x+y<1][i]
	T, save_interval, param_name, pert = (1000.0, 10.0, "cp_experiment", "double")
	param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert = pert)
	return param
end

function is_valid(f, param)
	return (mean(f[:, 1]) + mean(f[:, 2]) ≈ param["ϕa"]) & (mean(f[:, 3]) ≈ param["ϕp"])
end
#

## proccessing funcitons
function wave_speed(param::Dict{String, Any}, f::Matrix{Float64})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	nrho = sum(f; dims = 2)[:, 1]
	new_peak = argmax(nrho)
	old_peak = new_peak
	t = 0.0
	x_diff = 0
	#get first jump
	while new_peak == old_peak
		t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
		nrho = sum(f; dims = 2)[:, 1]
		new_peak = argmax(nrho)
	end
	# measure over 2 jumps
	t = 0.0
	old_peak = new_peak
	while x_diff<2
		t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR);
		nrho = sum(f; dims = 2)[:, 1]
		new_peak = argmax(nrho)
		x_diff = min(abs((new_peak-old_peak + Nx) % Nx), abs((old_peak-new_peak + Nx) % Nx))
	end
	c = Δx*x_diff/t

	t1, f1 = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR)
	t2, f2 = time_step!(t1, f1; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR)
	df = norm(f2-f)/sqrt(Nx)/(t2-t)
	d2f = norm((f2-f1)/(t2-t1)-(f1-f)/(t1-t))/sqrt(Nx)/(t2-t)/2

	filename = speed_save_name(param)
	data = Dict{String, Any}()
	@pack! data = c, df, d2f, f
	safesave(filename, data)

	return c, df, d2f, f
end

function n_step!(t, f, n; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR)
	for _ in 1:n
		t, f = time_step!(t, f; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR)
	end
	return t, f
end

function easy_speed(param::Dict{String, Any}, f::Matrix{Float64})
	@unpack DT, v0, DR, Δx, Lx, Nx, δt = param
	t1, f1 = n_step!(t, f, 10; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR)
	df = norm(f1-f)/sqrt(Nx)/(t1-t)
	fx = norm((f - circshift(f, (1, 0))))/sqrt(Nx)/Δx
	return df/fx
end

function f_dot(param::Dict{String, Any}, f::Matrix{Float64})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
	t = 0.0
	t1, f1 = n_step!(t, f, 10; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR)
	t2, f2 = n_step!(t1, f1, 10; δt = δt, Nx = Nx, Lx = Lx, DT = DT, v0 = v0, DR = DR)

	df1 = norm(f1-f)/sqrt(Nx)/(t1-t)
	df2 = norm(f2-f1)/sqrt(Nx)/(t2-t1)

	ddf = (df2 - df1)/(t2-t)/2

	fx1 = norm((f1 - circshift(f1, (1, 0))))/sqrt(Nx)/Δx
	fx2 = norm((f2 - circshift(f2, (1, 0))))/sqrt(Nx)/Δx

	c1 = df1/fx1
	c2 = df2/fx2

	dc = (c2 - c1)/(t2-t)/2

	normf = sqrt(sum((f[:, 1] .- ϕa/2) .^ 2 + (f[:, 2] .- ϕa/2) .^ 2 + (f[:, 3] .- ϕp) .^ 2)/Nx)

	return [normf, c1, dc]
end

function comp_wave_speed(param::Dict{String, Any})
	loaded, f, t = load_last_pde(param)
	if loaded
		return wave_speed(param, f)
	else
		return 0, 0, 1000.0, f
	end
end

function load_wave_speed(param::Dict{String, Any})
	try
		filename = speed_save_name(param)
		data::Dict{String, Any} = load(filename)
		@unpack c, df, d2f, f = data
		return c, df, d2f, f
	catch
		@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T, name, Nx, save_interval, save_on, δt = param
		println("load failed for [DT,v0,DR,Δx,Lx,ϕa,ϕp]=$([DT,v0,DR,Δx,Lx,ϕa,ϕp])")
		return 0.0, 0.0, 1000.0, initiate_uniform_pde(ϕa, ϕp, Nx);
	end
end

function t_dff(ts, ft; N = 100, gap = 1)
	ts = ts[gap:gap:end]
	Nt = length(ts)

	f_dt = zeros(Nt)

	for i in 2:Nt
		local fdiff, tdiff
		tdiff = ts[i] - ts[i-1]
		fdiff = ft[i*gap, :, :] - ft[(i-1)*gap, :, :]
		f_dt[i] = norm(fdiff / tdiff)/sqrt(N)
	end
	return ts, f_dt
end
#

## ode fns

function stretch_param(param, LLx)
	global param
	param["Δx"] = LLx/param["Nx"]
	param["Lx"] = LLx
	return param
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

function DD(ρ; logtol = 1e-10)
	α::Float64 = π/2 - 1;
	if ρ ≤ 0.0
		ρ = logtol
	elseif ρ>1.0
		ρ = 1.0
	end

	return -(α*(2*α-1)/(2*α+1)*ρ - α)+(α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ + 1)
end

function DDp(ρ; logtol = 1e-10)
	α::Float64 = π/2 - 1;
	if ρ ≤ 0.0
		ρ = logtol
	elseif ρ>1.0
		ρ = 1.0
	end

	return -(α*(2*α-1)/(2*α+1))+(2*α*(2*α-1)/(2*α+1)*ρ - α)
end

function self_diff(ρ; logtol = 1e-10)
	α = π/2 - 1;
	if ρ ≤ 0.0
		ρ = logtol
	elseif ρ>1.0
		ρ = 1.0
	end
	return (1-ρ) .* (α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ + 1)
end

function self_diff_prime(ρ; logtol = 1e-10)
	α = π/2 - 1;
	if ρ ≤ 0.0
		ρ = logtol
	elseif ρ>1.0
		ρ = 1.0
	end
	return - (α*(2*α-1)/(2*α+1)*ρ .^ 2 - α*ρ .+ 1) + (-ρ .+ 1)*(2*α*(2*α-1)/(2*α+1)*ρ - α);
end

ss(x) = DD(x) - 1
ds(x) = self_diff(x)
dsp(x) = self_diff_prime(x)
ip(f, i) = circshift(f, (-i, 0))


function dspp(ρ)
	α = π/2 - 1;
	return -2*(2*α*(2*α-1)/(2*α+1)*ρ - α) + (-ρ .+ 1)*(2*α*(2*α-1)/(2*α+1));
end
ssp(x) = DDp(x)


function get_f(u)
	Nx = (length(u)-1)÷3
	ρ = u[(0*Nx+1):1:(1*Nx)]
	ρa = u[(1*Nx+1):1:(2*Nx)]
	m = u[(2*Nx+1):1:(3*Nx)]
	c = u[end]
	f = zeros(Nx, 3)
	f[:, 2] = ρa/2 + m/2
	f[:, 1] = ρa/2 - m/2
	f[:, 3] = ρ - ρa
	return f
end

function get_u(f, c)
	# flatten initial guess
	ρ = sum(f; dims = 2)
	ρa = f[:, 2] + f[:, 1]
	m = f[:, 2] - f[:, 1]
	C = [c]

	return vcat(ρ, ρa, m, C)
end

function steady_save_name(param::Dict{String, Any})
	@unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, name, save_interval = param
	return datadir("pm_pdes_pro", "steady_state", "$(name)/[DT,v0,DR,Δx,Lx,ϕa,ϕp]=$([DT,v0,DR,Δx,Lx,ϕa,ϕp]).jld2")
end



# find outer soluton closest to a target global densities

function find_outer_sol(target_ϕa::Float64, target_ϕp::Float64, outer_type::Int64)
	# Get all .jld2 files in current directory

	if outer_type==0
		targetdir = datadir("pm_pdes_pro", "outer_sol", "outer_sol_0")
		pattern = r"\[DT,v0,DR,Nx,ϕa,ϕp\]=\[([^,]+),([^,]+),([^,]+),([^,]+),([^,]+),([^\]]+)\]\.jld2"

	elseif outer_type==1
		targetdir = datadir("pm_pdes_pro", "outer_sol", "outer_sol_1")
		pattern = r"\[DT,v0,DR,Nx,ϕ,γ\]=\[([^,]+),([^,]+),([^,]+),([^,]+),([^,]+),([^\]]+)\]\.jld2"

	elseif outer_type==2
		targetdir = datadir("pm_pdes_pro", "outer_sol", "outer_sol_2")
		pattern = r"\[DT,v0,DR,Nx,ind,γ\]=\[([^,]+),([^,]+),([^,]+),([^,]+),([^,]+),([^\]]+)\]\.jld2"

	else
		targetdir = datadir("pm_pdes_pro", "outer_sol");
	end

	files = filter(f -> endswith(f, ".jld2"), readdir(targetdir, join = true))
	println("Number of files in folder: ", length(files))

	# Initialize variables to track closest match
	closest_file = ""
	min_distance = Inf


	# Function to calculate Euclidean distance
	function calculate_distance(ϕa1, ϕp1, ϕa2, ϕp2)
		return sqrt((ϕa1 - ϕa2)^2 + (ϕp1 - ϕp2)^2)
	end

	# Iterate over files
	for file in files
		m = match(pattern, file)
		if !isnothing(m)
			# Extract parameters
			DT = parse(Float64, m[1])
			v0 = parse(Float64, m[2])
			DR = parse(Float64, m[3])
			Nx = parse(Float64, m[4])

			# Construct full filepath
			filename = joinpath(targetdir, file)
			data = load(filename)
			@unpack f = data

			ϕp = sum(f)/Nx-sum(f[:, 1:2])/Nx
			ϕa = sum(f[:, 1:2])/Nx


			distance = calculate_distance(ϕa, ϕp, target_ϕa, target_ϕp)

			# Update closest file if distance is smaller
			if distance < min_distance
				min_distance = distance
				closest_file = file
			end
		end
	end

	# Return results
	if closest_file != ""
		return (closest_file, min_distance)
	else
		return nothing
	end
end

#

## outer fns
function chop_f(f)
	Nx, _ = size(f)
	g = circshift(f, (Nx÷2, 0))
	rho = sum(g; dims = 2)[:, 1]
	I = argmin(rho):1:argmax(g[:, 1]+g[:, 2])
	return [g[i, j] for i in I, j in 1:3]
end

function resize_f(f, LLx, NNx)
	Nx, _ = size(f)
	Δx = 1/Nx
	ΔΔx = 1/NNx
	function fx(x, k)
		i::Int64 = x÷Δx
		j = i+1
		if i == 0
			i = 1
		end
		if j == Nx+1
			j = Nx
		end
		return ((Δx*(i+1) - x)*f[i, k] + (x-Δx*i)*f[j, k])/(Δx)
	end
	g = [fx(x, k) for x in ΔΔx:ΔΔx:1, k in 1:3]

	Δx = LLx/NNx
	Lx = LLx
	Nx = NNx
	@pack! param = Nx, Lx, Δx
	return param, circshift(g, (NNx÷2, 0))
end


function get_out_f(u, param; Lx = Inf)
	@unpack DT, v0, DR, Δx, Nx, ϕa, ϕp = param
	Nx = (length(u)-1)÷2
	ρ = u[(0*Nx+1):1:(1*Nx)]
	ρa = u[(1*Nx+1):1:(2*Nx)]

	# first order magnetisation
	m = -(1/Lx)*(v0/2/sqrt(DT*DR)/(2*Δx))*(circshift(ρa .* ds.(ρ), -1) - circshift(ρa .* ds.(ρ), 1))
	f = zeros(Nx, 3)
	f[:, 2] = ρa/2 + m/2
	f[:, 1] = ρa/2 - m/2
	f[:, 3] = ρ - ρa
	return circshift(f, (Nx÷2, 0))
end

function get_out_u(f, c, Lx; shift = true)
	# flatten initial guess
	Nx, _ = size(f)
	if shift
		f = circshift(f, (Nx÷2, 0))
	end
	ρ = sum(f; dims = 2)
	ρa = f[:, 2] + f[:, 1]
	C = [c*Lx]
	return vcat(ρ, ρa, C)
end

function get_out_u_fromfull(f, c, Lx; shift = true)
	Nx, _ = size(f)
	if shift
		f = circshift(f, (Nx÷2, 0))
	end
	ρ = sum(f; dims = 2)
	ρa = f[:, 2] + f[:, 1]
	C = [c*Lx]
	return vcat(ρ, ρa, C)
end

function get_BC(γ, Pe; initial_Δ = 1e-5, max_iter = 20, tol = 1e-8)
	find_sol = false
	try
		global lower_limits, upper_limits
		find_sol, lower_limits, upper_limits = colapse_sol_interval(; Pe = Pe, γ = γ, initial_Δ = initial_Δ, max_iter = max_iter, tol = tol)
	catch
		find_sol = false
	end
	if find_sol
		ϕL, ϕR = lower_limits[1], upper_limits[1]
		ϕaL, ϕaR = 1-γ*(1-ϕL), 1-γ*(1-ϕR)
		return ϕaL, ϕaR, ϕL, ϕR
	else
		print("no solution: output γ_max")
		ϕL = ϕR = find_rho_limit(; Pe = Pe, initial_Δ = initial_Δ, γ_max = 100.0)
		ϕaL, ϕaR = 1-γ*(1-ϕL), 1-γ*(1-ϕR)
		return ϕaL, ϕaR, ϕL, ϕR
	end
end

function initial_guess(ps)
	ps = (DT, v0, DR, Δx, Nx, Lx, ϕa, ϕaL, ϕaR, ϕL, ϕR) = ps
	xs = collect(0:(1/(Nx-1)):1)
	ρ = xs*(ϕR - ϕL) .+ ϕL
	ρa = xs*(ϕaR-ϕaL) .+ ϕaL
	C = [1.5]
	return vcat(ρ, ρa, C)
end

function hhh0(x::ForwardDiff.Dual{T, V, N}, γ::ForwardDiff.Dual{T, V, N}; Pe = 10, atol = 1e-12, log_tol = 1e-15) where {T, V, N}
	x_val = x.value
	y_val = γ.value
	if x_val > 1.0-log_tol
		x_val = 1-log_tol
	end
	fx = hhh0(x_val, y_val; Pe = Pe, atol = atol)
	dfdx_val = dgg0(x_val, y_val; Pe = Pe)*R.(x_val; tol = atol)
	dfdy_val = gg0_γ(x_val, y_val; Pe = Pe)*R.(x_val; tol = atol) - ΦoR_γ(x_val, Pe, y_val)
	return ForwardDiff.Dual{T}(fx, (dfdx_val*ForwardDiff.partials(x) + dfdy_val*ForwardDiff.partials(γ)))
end

function hhh0(ρ::Float64, γ::Float64; Pe = 10, atol = 1e-12, log_tol = 1e-15)
	if ρ > 1.0-log_tol
		ρ = 1-log_tol
	end
	return gg0(ρ, γ; Pe = Pe)*R.(ρ; tol = atol)-ΦoR.(ρ; Pe = Pe, γ = γ, tol = atol)
end

function DS(ρ; logtol = 1e-10)
	α = π/2 - 1;
	return (1-ρ) .* (α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ + 1)
end

function gg0(ρ, γ; Pe = 10, log_tol = 1e-15)
	if ρ > 1.0-log_tol
		ρ = 1-log_tol
	end
	return Pe .* (-γ*(-ρ .+ 1) .+ 1) .* DS.(ρ) - 2*log.(- ρ .+ 1) ./ Pe
end
dgg0(ρ, γ; Pe = 10) = Pe*(γ)*DS.(ρ) + Pe*(-γ*(- ρ .+ 1) .+ 1)*dsp.(ρ) + 2/(- ρ .+ 1)/Pe
gg0_γ(ρ, γ; Pe = 10) = Pe*(-(1 - ρ))*DS(ρ)
# gg0(x::ForwardDiff.Dual{T,V,N};Pe = 10, γ::ForwardDiff.Dual{T,V,N} = 1, atol = 1 ) where {T,V,N} = ForwardDiff.Dual{T}( g0(x.value;Pe = Pe, γ = γ), dg0(x.value;Pe = Pe, γ = γ)*ForwardDiff.partials(x) )

dΦ_dρ_γ(ρ, Pe, γ) = -Pe/(1-ρ)^2
ΦoR_γ(ρ, Pe, γ) = -Pe*(-ρ .+ 1) .^ (-1) .+ Pe



function increase_ind(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind = ps
	i = ind
	u[i+1] = u[i]
	i = Nx + i
	u[i+1] = u[i]
	ps = DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, (ind+1)
	return u, ps
end

function decrease_ind(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind = ps
	i = ind
	u[i] = u[i+1]
	i = Nx + i
	u[i] = u[i+1]
	ps = DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind-1
	return u, ps
end

function extract_ind(f)
	Nx, _ = size(f)
	m = f[:, 2] - f[:, 1];
	return argmax(m) + Nx÷2
end
#

##############################
## FULL TRAVELLING SOLUTION
##############################


function rho_eq(ρ, ρa, m, c, i; param = p)
	DT, v0, DR, Δx, Nx, Lx = param
	iplus = (i + Nx)%Nx + 1
	iminus = (i + Nx - 2)%Nx + 1
	return -c*ρ[i] + (DT/2/Δx)*(ρ[iplus]-ρ[iminus]) - v0*(1-ρ[i])*m[i]
end

function act_eq(ρ, ρa, m, c, i; param = p)
	DT, v0, DR, Δx, Nx, Lx = param
	iplus = (i + Nx)%Nx + 1
	iminus = (i + Nx - 2)%Nx + 1
	return -c*ρa[i] + (DT/2/Δx)*(ds(ρ[i])*(ρa[iplus]-ρa[iminus]) + ρa[i]*DD(ρ[i])*(ρ[iplus]-ρ[iminus])) - v0*(ρa[i]*ss(ρ[i])*m[i] + ds(ρ[i])*m[i])
end

function mag_eq(ρ, ρa, m, c, i; param = p)
	DT, v0, DR, Δx, Nx, Lx = param
	iplus = (i + Nx)%Nx + 1
	iminus = (i + Nx - 2)%Nx + 1
	return -c*m[i] + (DT/2/Δx)*(ds(ρ[i])*(m[iplus]-m[iminus]) + m[i]*DD(ρ[i])*(ρ[iplus]-ρ[iminus])) - v0*(ss(ρ[i])*m[i]^2 + ds(ρ[i])*ρa[i])
end

function ff(F, u, p)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp = p
	ϕ = ϕa + ϕp
	global ρ, ρa, m, c
	ρ  = u[(0*Nx+1):1:(1*Nx)]
	ρa = u[(1*Nx+1):1:(2*Nx)]
	m   = u[(2*Nx+1):1:(3*Nx)]
	c   = u[end]

	for i in 1:(Nx-1)
		F[i]        = rho_eq(ρ, ρa, m, c, i; param = p) - rho_eq(ρ, ρa, m, c, Nx; param = p)
		F[Nx-1+i]   = act_eq(ρ, ρa, m, c, i; param = p) - act_eq(ρ, ρa, m, c, Nx; param = p)
		F[2*Nx-2+i] = mag_eq(ρ, ρa, m, c, i; param = p) - mag_eq(ρ, ρa, m, c, Nx; param = p) - 2*DR*Δx*sum(m[1:i])
	end

	F[3*Nx-2] = sum(m)/Nx
	F[3*Nx-1] = sum(ρ)/Nx - ϕ
	F[3*Nx]   = sum(ρa)/Nx - ϕa
	F[3*Nx+1] = ρ[end] - ϕ
end

function return_f(u; p = ps)
	F   = deepcopy(u)
	ρ  = u[(0*Nx+1):1:(1*Nx)]
	ρa = u[(1*Nx+1):1:(2*Nx)]
	m   = u[(2*Nx+1):1:(3*Nx)]
	c   = u[end]

	for i in 1:(Nx-1)
		F[i]        = rho_eq(ρ, ρa, m, c, i; param = p) - rho_eq(ρ, ρa, m, c, Nx; param = p)
		F[Nx-1+i]   = act_eq(ρ, ρa, m, c, i; param = p) - act_eq(ρ, ρa, m, c, Nx; param = p)
		F[2*Nx-2+i] = mag_eq(ρ, ρa, m, c, i; param = p) - mag_eq(ρ, ρa, m, c, Nx; param = p) - 2*DR*Δx*sum(m[1:i])
	end

	F[3*Nx-2] = sum(m)/Nx
	F[3*Nx-1] = sum(ρ)/Nx - ϕ
	F[3*Nx]   = sum(ρa)/Nx - ϕa
	F[3*Nx+1] = ρ[end] - ϕ
	return F
end

function load_full(Lx, Nx, ϕa, ϕp, v0)
	param, ps = get_param_full(Lx, Nx, ϕa, ϕp, v0)
	filename = save_name_full(param)
	data = load(filename)
	println(filename)
	@unpack f, c = data
	u = get_u(f, c)
	return f, u, c
end

function solve_full(Lx, Nx, ϕa, ϕp, v0, u0; tol = 1e-8, maxiters = 100)
	param, p = get_param_full(Lx, Nx, ϕa, ϕp, v0)
	# Set problem
	prob = NonlinearProblem(ff, u0, p; abstol = tol, reltol = tol, maxiters = maxiters);
	# sol  = solve(prob, callback = cb)
	sol = solve(prob)
	#
	# Save
	u        = sol.u
	f        = get_f(u)
	c        = u[end]
	filename = save_name_full(param)
	data     = Dict("f" => f, "c" => c)
	safesave(filename, data)
	#
	return f, u, c
end

function get_param_full(Lx, Nx, ϕa, ϕp, v0)
	param = get_stretch_param(Lx)
	name = "sol_full"
	Δx = Lx/Nx
	DT = DR = 1.0
	@pack! param = Lx, Nx, ϕa, ϕp, v0, name, Δx
	ps = DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp
	return param, ps
end

function save_name_full(param)
	@unpack DT, v0, DR, Nx, Lx, name, ϕa, ϕp = param
	return datadir("pm_pdes_pro", "full_sol", "$(name)/[DT,v0,DR,Lx,Nx,ϕa,ϕp]=$([DT,v0,DR,Lx,Nx,ϕa,ϕp]).jld2")
end

function check_u_full(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp = ps
	param, _ = get_param_full(Lx, Nx, ϕa, ϕp, v0)
	F = zeros(3*Nx+1)
	ff(F, u, ps)
	f = get_f(u)
	avmag = Lx*sum(f[:, 2]-f[:, 1])/Nx
	return maximum(abs.(F)), argmax(abs.(F)), avmag, u[end]
end

function align_f(f, c)
	ρ = sum(f; dims = 2)[:, 1]
	ϕ = mean(ρ)
	i = argmax(ρ)
	ρ = circshift(ρ, -i)
	k = argmin(ρ)
	j = argmin(abs.(ρ[1:k] .- ϕ))
	f = circshift(f, (-i-j, 0))
	u = get_u(f, c)
	return f, u, c
end

#############################
# OUTER REGION EQUATIONS
#############################

function rho_out_eq(ρ, ρa, c, i; p = p)
	DT, v0, DR, _, Nx, _, ϕa, ϕp = p
	Δx = 1/Nx
	Lx = 1
	iplus = i+1
	ρ1 = ρ[i]
	ρ2 = ρ[i+1]
	ρ3 = 2*(ρ1*ρ2)/(ρ1+ρ2)
	ρa1 = ρa[i]
	ρa2 = ρa[i+1]
	ρa3 = 2*(ρa1*ρa2)/(ρa1+ρa2)
	return -c*(ρ2+ρ1)/2 + (DT/Δx)*(ρ2-ρ1) + (v0^2/DR/Δx/2)*(ds(ρ3)*(1-ρ3)*(ρa2-ρa1) + ρa3*dsp(ρ3)*(1-ρ3)*(ρ2-ρ1))
end

function act_out_eq(ρ, ρa, c, i; p = p)
	DT, v0, DR, _, Nx, _, ϕa, ϕp = p
	Δx = 1/Nx
	Lx = 1
	iplus = i + 1
	ρ1 = ρ[i]
	ρ2 = ρ[i+1]
	ρ3 = 2*(ρ1*ρ2)/(ρ1+ρ2)
	ρa1 = ρa[i]
	ρa2 = ρa[i+1]
	ρa3 = 2*(ρa1*ρa2)/(ρa1+ρa2)
	return -c*(ρa2+ρa1)/2 + (DT/Δx)*(ds(ρ3)*(ρa2-ρa1) + ρa3*DD(ρ3)*(ρ2-ρ1)) + (v0^2/DR/Δx/2)*(ds(ρ3)*(ρa3*ss(ρ3)+ds(ρ3))*(ρa2-ρa1) + ρa3*dsp(ρ3)*(ρa3*ss(ρ3)+ds(ρ3))*(ρ2-ρ1))
end

##########################
## TRAVELLING SOLUTION, OUTER PROBLEM WITH NO INTERFACES, OUTER0
##########################


function ff_outer0(F, u, p; atol = 1e-12)
	DT, v0, DR, _, Nx, _, ϕa, ϕp = p
	global ρ, ρa, m
	local c
	ρ  = u[(0*Nx+1):1:(1*Nx)]
	ρa = u[(1*Nx+1):1:(2*Nx)]
	c   = u[end]

	for i in 2:(Nx-1)
		F[i-1]    = rho_out_eq(ρ, ρa, c, i; p = p) - rho_out_eq(ρ, ρa, c, 1; p = p)
		F[Nx-3+i] = act_out_eq(ρ, ρa, c, i; p = p) - act_out_eq(ρ, ρa, c, 1; p = p)
	end

	Δx = 1/Nx
	ρ1 = ρ[end]
	ρ2 = ρ[1]
	ρ3 = 2*(ρ1*ρ2)/(ρ1+ρ2)
	ρa1 = ρa[end]
	ρa2 = ρa[1]
	ρa3 = 2*(ρa1*ρa2)/(ρa1+ρa2)

	F[2*Nx-3] = -c*(ρ2+ρ1)/2 + (DT/Δx)*(ρ2-ρ1) + (v0^2/DR/Δx/2)*(ds(ρ3)*(1-ρ3)*(ρa2-ρa1) + ρa3*dsp(ρ3)*(1-ρ3)*(ρ2-ρ1)) - rho_out_eq(ρ, ρa, c, 1; p = p)
	F[2*Nx-2] = -c*(ρa2+ρa1)/2 + (DT/Δx)*(ds(ρ3)*(ρa2-ρa1) + ρa3*DD(ρ3)*(ρ2-ρ1)) + (v0^2/DR/Δx/2)*(ds(ρ3)*(ρa3*ss(ρ3)+ds(ρ3))*(ρa2-ρa1) + ρa3*dsp(ρ3)*(ρa3*ss(ρ3)+ds(ρ3))*(ρ2-ρ1)) - act_out_eq(ρ, ρa, c, 1; p = p)
	F[2*Nx-1] = ρ[Nx÷2] - ϕa - ϕp
	F[2*Nx]   = sum(ρ-ρa)/Nx - ϕp
	F[2*Nx+1] = sum(ρa)/Nx - ϕa
end

function load_outer0(Lx, Nx, ϕa, ϕp, v0)
	param, ps = get_param_outer0(Lx, Nx, ϕa, ϕp, v0)
	filename = save_name_outer0(param)
	data = load(filename)
	println(filename)
	@unpack f, c = data
	u = get_out_u(f, c, Lx)
	return f, u, c
end


function solve_outer0(Lx, Nx, ϕa, ϕp, v0, u0; tol = 1e-8, maxiters = 100)
	param, p = get_param_outer0(Lx, Nx, ϕa, ϕp, v0)
	# Set problem
	prob = NonlinearProblem(ff_outer0, u0, p; abstol = tol, reltol = tol, maxiters = maxiters);
	sol  = solve(prob)
	#
	# Save
	u        = sol.u
	f        = get_out_f(u, param)
	c        = u[end]/Lx
	filename = save_name_outer0(param)
	data     = Dict("f" => f, "c" => c)
	safesave(filename, data)
	#
	return f, u, c
end

function get_param_outer0(Lx, Nx, ϕa, ϕp, v0)
	param = get_stretch_param(Lx)
	name = "outer_sol_0"
	Δx = Lx/Nx
	DT = DR = 1.0
	@pack! param = Lx, Nx, ϕa, ϕp, v0, name, Δx
	ps = (DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp)
	return param, ps
end

function save_name_outer0(param)
	@unpack DT, v0, DR, Nx, Lx, name, ϕa, ϕp = param
	return datadir("pm_pdes_pro", "outer_sol", "$(name)/[DT,v0,DR,Nx,ϕa,ϕp]=$([DT,v0,DR,Nx,ϕa,ϕp]).jld2")
end



function check_u_outer0(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp = ps
	param, _ = get_param_outer0(Lx, Nx, ϕa, ϕp, v0)
	F = zeros(2*Nx+1)
	ff_outer0(F, u, ps)
	f = get_out_f(u, param)
	avmag = Lx*sum(f[:, 2]-f[:, 1])/Nx
	return maximum(abs.(F)), avmag, u[end]
end
#

##########################
## TRAVELLING SOLUTION, OUTER PROBLEM WITH ONE INTERFACE, OUTER1
##########################
function ff_outer1(F, u, p; atol = 1e-12)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ϕ, γ, ϕaL, ϕaR, ϕL, ϕR = p
	global ρ, ρa, m, ep
	local c
	ρ  = u[(0*Nx+1):1:(1*Nx)]
	ρa = u[(1*Nx+1):1:(2*Nx)]
	c   = u[end]

	for i in 2:(Nx-1)
		F[i-1]    = rho_out_eq(ρ, ρa, c, i; p = p) - rho_out_eq(ρ, ρa, c, 1; p = p)
		F[Nx-3+i] = act_out_eq(ρ, ρa, c, i; p = p) - act_out_eq(ρ, ρa, c, 1; p = p)
	end

	F[2*Nx-3] = ρa[1] - ϕaL
	F[2*Nx-2] = ρa[end] - ϕaR
	F[2*Nx-1] = ρ[1] - ϕL
	F[2*Nx]   = ρ[end] - ϕR
	F[2*Nx+1] = sum(ρ)/Nx - ϕ
end

function load_outer1(Lx, Nx, ϕa, ϕp, v0, ϕ, γ)
	param, ps = get_param_outer1(Lx, Nx, ϕa, ϕp, v0, ϕ, γ)
	filename = save_name_outer1(param, ϕ, γ)
	data = load(filename)
	println(filename)
	@unpack f, c = data
	u = get_out_u(f, c, Lx)
	return f, u, c
end

function solve_outer1(Lx, Nx, ϕa, ϕp, v0, ϕ, γ, u0; tol = 1e-8, maxiters = 100)
	param, p = get_param_outer1(Lx, Nx, ϕa, ϕp, v0, ϕ, γ)
	# Set problem
	prob = NonlinearProblem(ff_outer1, u0, p; abstol = tol, reltol = tol, maxiters = maxiters);
	sol  = solve(prob)
	#
	# Save
	u        = sol.u
	f        = get_out_f(u, param)
	c        = u[end]/Lx
	filename = save_name_outer1(param, ϕ, γ)
	data     = Dict("f" => f, "c" => c)
	safesave(filename, data)
	#
	return f, u, c
end

function get_param_outer1(Lx, Nx, ϕa, ϕp, v0, ϕ, γ)
	param = get_stretch_param(Lx)
	name = "outer_sol_1"
	Δx = Lx/Nx
	DT = DR = 1.0
	ϕaL, ϕaR, ϕL, ϕR = get_BC(γ, v0; initial_Δ = 1e-6, tol = 1e-8)
	@pack! param = Lx, Nx, ϕa, ϕp, v0, name, Δx, ϕ, γ
	ps = DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ϕ, γ, ϕaL, ϕaR, ϕL, ϕR
	return param, ps
end

function save_name_outer1(param, ϕ, γ)
	@unpack DT, v0, DR, Nx, Lx, name, ϕa, ϕp = param
	return datadir("pm_pdes_pro", "outer_sol", "$(name)/[DT,v0,DR,Nx,ϕ,γ]=$([DT,v0,DR,Nx,ϕ,γ]).jld2")
end

function check_u_outer1(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ϕ, γ, ϕaL, ϕaR, ϕL, ϕR = ps
	param, _ = get_param_outer1(Lx, Nx, ϕa, ϕp, v0, ϕ, γ)
	F = zeros(2*Nx+1)
	ff_outer1(F, u, ps)
	f = get_out_f(u, param)
	avmag = Lx*sum(f[:, 2]-f[:, 1])/Nx
	return maximum(abs.(F)), argmax(abs.(F)), avmag, u[end]/Lx
end
#


##########################
## TRAVELLING SOLUTION, OUTER PROBLEM WITH TWO INTERFACES, OUTER2
##########################

function ff_outer2(F, u, p; atol = 1e-12)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind, γ, ϕaL, ϕaR, ϕL, ϕR = p
	global ρ, ρa, m, ep
	local c
	ρ  = u[(0*Nx+1):1:(1*Nx)]
	ρa = u[(1*Nx+1):1:(2*Nx)]
	c   = u[end]

	for i in 2:(ind-1)
		F[i-1]    = rho_out_eq(ρ, ρa, c, i; p = p) - rho_out_eq(ρ, ρa, c, 1; p = p)
		F[Nx-3+i] = act_out_eq(ρ, ρa, c, i; p = p) - act_out_eq(ρ, ρa, c, 1; p = p)
	end

	for i in (ind+1):1:(Nx-1)
		F[i-1]    = rho_out_eq(ρ, ρa, c, i; p = p) - rho_out_eq(ρ, ρa, c, 1; p = p)
		F[Nx-3+i] = act_out_eq(ρ, ρa, c, i; p = p) - act_out_eq(ρ, ρa, c, 1; p = p)
	end

	ind1        = ind + 1
	γ_end      = (-ρa[ind] .+ 1) ./ (-ρ[ind] .+ 1)
	γ_srt      = (-ρa[ind1] .+ 1) ./ (-ρ[ind1] .+ 1)
	F[ind-1]    = gg0.(ρ[ind], γ_end; Pe = v0) - gg0.(ρ[ind1], γ_srt; Pe = v0)
	F[Nx-3+ind] = (hhh0.(ρ[ind], γ_end; Pe = v0, atol = atol) - hhh0.(ρ[ind1], γ_srt; Pe = v0, atol = atol)) ./ 10
	F[2*Nx+1]   = γ_end .- γ_srt

	F[2*Nx-3] = ρa[1] - ϕaL
	F[2*Nx-2] = ρa[end] - ϕaR
	F[2*Nx-1] = ρ[1] - ϕL
	F[2*Nx]   = ρ[end] - ϕR
end

function load_outer2(Lx, Nx, ϕa, ϕp, v0, ind, γ)
	param, ps = get_param_outer2(Lx, Nx, ϕa, ϕp, v0, ind, γ)
	filename = save_name_outer2(param, ind, γ)
	data = load(filename)
	println(filename)
	@unpack f, c = data
	u = get_out_u(f, c, Lx)
	return f, u, c
end

function solve_outer2(Lx, Nx, ϕa, ϕp, v0, ind, γ, u0; tol = 1e-8, maxiters = 100)
	param, p = get_param_outer2(Lx, Nx, ϕa, ϕp, v0, ind, γ)
	# Set problem
	prob = NonlinearProblem(ff_outer2, u0, p; abstol = tol, reltol = tol, maxiters = maxiters);
	sol  = solve(prob)
	#
	# Save
	u        = sol.u
	f        = get_out_f(u, param)
	c        = u[end]/Lx
	filename = save_name_outer2(param, ind, γ)
	data     = Dict("f" => f, "c" => c)
	safesave(filename, data)
	#
	return f, u, c
end

function get_param_outer2(Lx, Nx, ϕa, ϕp, v0, ind, γ)
	param = get_stretch_param(Lx)
	name = "outer_sol_2"
	Δx = Lx/Nx
	DT = DR = 1.0
	ϕaL, ϕaR, ϕL, ϕR = get_BC(γ, v0; initial_Δ = 1e-6, tol = 1e-8)
	@pack! param = Lx, Nx, ϕa, ϕp, v0, name, Δx, ind, γ
	ps = DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind, γ, ϕaL, ϕaR, ϕL, ϕR
	return param, ps
end

function save_name_outer2(param, ind, γ)
	@unpack DT, v0, DR, Nx, Lx, name, ϕa, ϕp = param
	return datadir("pm_pdes_pro", "outer_sol", "$(name)/[DT,v0,DR,Nx,ind,γ]=$([DT,v0,DR,Nx,ind,γ]).jld2")
end

function check_u_outer2(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind, γ, ϕaL, ϕaR, ϕL, ϕR = ps
	param, _ = get_param_outer2(Lx, Nx, ϕa, ϕp, v0, ind, γ)
	F = zeros(2*Nx+1)
	ff_outer2(F, u, ps)
	f = get_out_f(u, param)
	avmag = Lx*sum(f[:, 2]-f[:, 1])/Nx
	return maximum(abs.(F)), argmax(abs.(F)), avmag, u[end]/Lx
end

function increase_ind_outer2(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind, γ, ϕaL, ϕaR, ϕL, ϕR = ps
	i = ind
	u[i+1] = u[i]
	i = Nx + i
	u[i+1] = u[i]
	ps = DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, (ind+1), γ, ϕaL, ϕaR, ϕL, ϕR
	return u, ps
end

function decrease_ind_outer2(u, ps)
	DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, ind, γ, ϕaL, ϕaR, ϕL, ϕR = ps
	i = ind
	u[i] = u[i+1]
	i = Nx + i
	u[i] = u[i+1]
	ps = DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp, (ind-1), γ, ϕaL, ϕaR, ϕL, ϕR
	return u, ps
end

function increase_inds(u, ps, n)
	for i in 1:n
		u, ps = increase_ind_outer2(u, ps)
	end
	return u, ps
end

function decrease_inds(u, ps, n)
	for i in 1:n
		u, ps = decrease_ind_outer2(u, ps)
	end
	return u, ps
end

