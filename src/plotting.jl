# plotting functions 

using KernelDensity, KernelDensitySJ, Peaks, Statistics, ForwardDiff;
using PyPlot, LaTeXStrings


function plot_sim(sim_ts, η_saves, param)
	# Process data for visualization
	ϵ = 0.05  # Smoothing parameter for local averaging

	images = []
	times  = []
	fts    = []

	# Process simulation data
	@unpack N, Nx, N₁, N₂, Lx, Ly = param
	ft_sim = local_average_timeseries(η_saves, ϵ, N, N₁, N₂)
	t_sim_rgb_image = rho_to_rgb(ft_sim)

	# Center simulation pattern
	pk = find_xpeak_ft(sim_ts, ft_sim; time_length = 0.1)
	centre = N₁ ÷ 2 + 1; # place peak at center of domain
	ft_sim = circshift(ft_sim, (0, -pk + centre, 0))
	t_sim_rgb_image = circshift(t_sim_rgb_image, (pk - centre, 0, 0))

	push!(images, t_sim_rgb_image)
	push!(times, sim_ts)
	push!(fts, ft_sim)


	# Define plot limits and layout parameters
	rhomin, rhomax = 0.0, 1.0
	mag_lim = 0.8
	font = 12  # Define font size

	# Layout dimensions
	height_1 = 0.08
	width_1 = 0.4
	side_gap_1 = 0.1
	bottom_gap_1 = 0.65
	gap = 0.012

	height_2 = 0.175
	width_2 = height_2
	side_gap_2 = 0.1
	bottom_gap_2 = 0.4
	gap_2 = 0.06

	# Annotation positions
	t_stamp_x = 0.03
	t_stamp_y = 0.05

	# Colorbar positions
	cbar_y_top = bottom_gap_1 + 3*height_1 + 2*gap
	cbar_width = 0.1
	cbar_y_bot = bottom_gap_1 + gap + height_1
	cbar_x = 2*side_gap_1 + width_1

	# Initialize figure
	rc("text", usetex = true)
	fig = plt.figure(figsize = (10, 10))

	# Create base axes for the figure
	ax = fig.add_axes([0, 0, 1, 1], visible = true)
	ax.xaxis.set_ticks([])
	ax.yaxis.set_ticks([])
	for spine in ["top", "right", "bottom", "left"]
		ax.spines[spine].set_visible(false)
	end

	@unpack T = param

	# Add image plots
	for (i, (rgb_image, ts)) in enumerate(zip(images, times))
		ax = fig.add_axes([side_gap_1+(i-1)*(side_gap_1+width_1), bottom_gap_1, width_1, height_1])
		t_end = ts[end]
		t_start = ts[1]
		ax.imshow(rgb_image; extent = [t_start, t_end, 0, Lx], interpolation = "bilinear")
		# ax.get_yaxis().set_ticks(0:1.0:Lx)
		# ax.set_yticklabels(["0","1","2"])
		# ax.get_xaxis().set_ticks([])
		# ax.axis([0,T, 0, 2])
		ax.set_aspect((T/2)*(height_1/width_1))
		ax.set_ylabel(L"x", fontsize = font, rotation = 90)
		# ax.get_xaxis().set_ticks(0:round(0.25*round(T); digits = 0):round(T;digits = 0 ))
		ax.set_xlabel(L"t", fontsize = font)
		ax.tick_params(labelbottom = true, direction = "in")
	end

	# Add matrix plots
	for (i, (ts, ft)) in enumerate(zip(times, fts))
		global im1, im2
		# Density plot
		ax = fig.add_axes([side_gap_1+(i-1)*(side_gap_1+width_1), bottom_gap_1+2*height_1+2*gap, width_1, height_1])
		t_end = ts[end]
		t_start = ts[1]
		_, Nx, _ = size(ft)

		colmap = PyPlot.plt.cm.viridis
		norm1 = matplotlib.colors.Normalize(vmin = rhomin, vmax = rhomax)
		im1 = ax.matshow((ft[:, Nx:-1:1, 1]+ft[:, Nx:-1:1, 2]+ft[:, Nx:-1:1, 3])';
			norm = norm1, cmap = colmap, extent = [t_start, t_end, 0, Lx])

		# ax.xaxis.set_ticks([])
		ax.xaxis.tick_bottom()
		# ax.get_yaxis().set_ticks(0:1.0:Lx)
		# ax.set_yticklabels(["0","1","2"])
		# ax.axis([0,T, 0, 2])
		ax.set_aspect((T/2)*(height_1/width_1))
		ax.set_ylabel(L"x", fontsize = font, rotation = 90)
		# ax.get_xaxis().set_ticks(0:round(0.25*round(T); digits = 0):round(T;digits = 0 ))
		ax.tick_params(labelbottom = false, direction = "in")

		# Magnetization plot
		ax = fig.add_axes([side_gap_1+(i-1)*(side_gap_1+width_1), bottom_gap_1+height_1+gap, width_1, height_1])
		colmap = PyPlot.plt.cm.PRGn
		norm1 = matplotlib.colors.Normalize(vmin = -mag_lim, vmax = mag_lim)
		im2 = ax.matshow((ft[:, Nx:-1:1, 2]-ft[:, Nx:-1:1, 1])';
			norm = norm1, cmap = colmap, extent = [t_start, t_end, 0, Lx])

		# ax.get_xaxis().set_ticks([])
		# ax.get_yaxis().set_ticks(0:1.0:Lx)
		# ax.set_yticklabels(["0","","2"])
		ax.xaxis.tick_bottom()
		# ax.axis([0,T, 0, 2])
		ax.set_aspect((T/2)*(height_1/width_1))
		ax.set_ylabel(L"x", fontsize = font, rotation = 90)
		# ax.get_xaxis().set_ticks(0:round(0.25*round(T); digits = 0):round(T;digits = 0 ))
		ax.tick_params(labelbottom = false, direction = "in")
	end

	# Add colorbars
	cbar_ax = fig.add_axes([cbar_x, bottom_gap_1, height_1, height_1])
	Δ = 0.001
	cbar_f = [x*(x+y≤1)*(i!=3)/2 + y*(x+y≤1)*(i==3) for x in Δ:Δ:1, y in Δ:Δ:1, i in 1:3]
	rgb_image = rho_to_rgb(cbar_f)

	ax = cbar_ax
	ax.imshow(rgb_image; extent = [0, 1, 0, 1])
	ax.spines["top"].set_visible(false)
	ax.spines["right"].set_visible(false)
	# ax.get_xaxis().set_ticks(0:0.5:1)
	# ax.get_yaxis().set_ticks(0:0.5:1)
	ax.set_xlabel(L"\rho_a", fontsize = font)
	ax.set_ylabel(L"\rho_0", fontsize = font, rotation = 90)
	ax.tick_params(direction = "in")

	# Add density colorbar
	rho_cbar_ax = fig.add_axes([cbar_x, cbar_y_bot+gap+height_1, 0.025, height_1])
	rho_cbar = fig.colorbar(im1, cax = rho_cbar_ax)
	rho_cbar_ax.set_ylabel(L"\rho", fontsize = font, rotation = 90)
	rho_cbar.set_ticks(rhomin:0.5:rhomax)
	rho_cbar_ax.yaxis.set_ticks_position("right")
	rho_cbar_ax.tick_params(direction = "in")

	# Add magnetization colorbar
	mag_cbar_ax = fig.add_axes([cbar_x, cbar_y_bot, 0.025, height_1])
	mag_cbar = fig.colorbar(im2, cax = mag_cbar_ax)
	mag_cbar.set_ticks((-mag_lim):0.8:mag_lim)
	mag_cbar_ax.tick_params(direction = "in")
	mag_cbar_ax.yaxis.set_ticks_position("right")
	mag_cbar_ax.set_ylabel(L"m", fontsize = font, rotation = 90)

	# Add final time plots
	# Total density (ρ)
	ax = fig.add_axes([side_gap_2, bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 2] + fts[1][end, :, 1] + fts[1][end, :, 3]
	# pde_mag = fts[2][end,:,2] + fts[2][end,:,1] + fts[2][end,:,3]

	@unpack N, Δx = param
	# ax.plot(Δx:Δx:Lx, pde_mag; color = "black", label = "PDE")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red", label = "Simulation")

	# ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks(rhomin:0.5:rhomax)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"\rho", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(rhomax-rhomin)))
	ax.axis([0, Lx, rhomin, rhomax])
	ax.tick_params(direction = "in")

	# Active density (ρₐ)
	ax = fig.add_axes([side_gap_2+gap_2+width_2, bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 2] + fts[1][end, :, 1]
	# pde_mag = fts[2][end,:,2] + fts[2][end,:,1]

	# ax.plot(Δx:Δx:Lx, pde_mag; color = "black")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red")

	# ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks(rhomin:0.5:rhomax)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"\rho_a", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(rhomax-rhomin)))
	ax.axis([0, Lx, rhomin, rhomax])
	ax.tick_params(direction = "in")

	# Passive density (ρ₀)
	ax = fig.add_axes([side_gap_2+2*(gap_2+width_2), bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 3]
	# pde_mag = fts[2][end,:,3] 

	# ax.plot(Δx:Δx:Lx, pde_mag; color = "black")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red")

	# ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks(rhomin:0.5:rhomax)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"\rho_0", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(rhomax-rhomin)))
	ax.axis([0, Lx, rhomin, rhomax])
	ax.tick_params(direction = "in")

	# Magnetization (m)
	ax = fig.add_axes([side_gap_2+3*(gap_2+width_2), bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 2] - fts[1][end, :, 1]
	# pde_mag = fts[2][end,:,2] - fts[2][end,:,1]

	# ax.plot(Δx:Δx:Lx, pde_mag; color = "black")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red")

	# ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks((-mag_lim):0.4:mag_lim)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"m", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(2*mag_lim)))
	ax.axis([0, Lx, -mag_lim, mag_lim])
	ax.tick_params(direction = "in")

	# Add time stamp and labels
	latex_annotation = latexstring("\$ t = $(round(times[1][end];digits = 1))\$")
	ax.annotate(latex_annotation, (t_stamp_x, bottom_gap_2+t_stamp_y),
		xycoords = "figure fraction", rotation = 90, fontsize = font)

	return fig;
end


function plot_pde(param, f)
	@unpack Δx, Lx, Nx = param

	fig, ax = subplots(1, 1, figsize = (10, 5))

	xs = Δx:Δx:Lx;

	ax.plot(xs, f[:, 1], label = L"\varrho_+");
	ax.plot(xs, f[:, 2], label = L"\varrho_-", linestyle = "--");
	ax.plot(xs, f[:, 3], label = L"\varrho_0");
	ax.plot(xs, f[:, 1]+f[:, 2]+f[:, 3], label = L"\varrho");
	ax.legend();

	return fig;
end

function plot_pde_outer(param, f)
	@unpack Δx, Lx, Nx = param

	fig, ax = subplots(1, 1, figsize = (10, 5))

	xs = Δx:Δx:Lx;

	ax.plot(xs, f[:, 1] + f[:, 2], label = L"\varrho_a");
	ax.plot(xs, f[:, 3], label = L"\varrho_0");
	ax.plot(xs, f[:, 1]+f[:, 2]+f[:, 3], label = L"\varrho");
	ax.legend();

	return fig;
end



"""
	plot_phase_reduced(fig, ax, Pe, font; Δϕ=0.01, plot_tieline = true, shading = true)

Plots the phase diagram with binodal and spinodal lines.

Arguments:
- fig: Figure object
- ax: Axis to plot on
- Pe: Peclet number
- font: Font size
- Δϕ: Phase space resolution (default 0.01)
- plot_tieline: Whether to plot the tie lines (default true)
- shading: Whether to shade the regions (default true)
"""
function plot_phase_reduced(fig, ax, Pe, font; Δϕ = 0.01, plot_tieline = true, shading = true)
	# Load binodal data
	filename = datadir("binodal", "Pe=$(Pe).jld2")
	data = wload(filename)
	@unpack Pe, γs, ϕ1s, ϕ2s = data

	# Set up axis labels and styling
	rc("text", usetex = true)
	ax.xaxis.set_tick_params(labelsize = font)
	ax.yaxis.set_tick_params(labelsize = font)
	ax.set_xlabel(L"\phi_a", fontsize = font*1.2)
	ax.set_ylabel(L"\phi_p", fontsize = font*1.2)
	ax.tick_params(labelbottom = true, direction = "in")

	# Plot spinodal lines
	ϕas_left, ϕas_right, ϕps, indl, indr = return_spin(; Pe = Pe, Δϕ = Δϕ)
	ax.plot(ϕas_left, ϕps, color = "blue", label = L"\mathrm{spinodal}", linestyle = "-")
	ax.plot(ϕas_right, ϕps, color = "blue", label = "_spinodal", linestyle = "-")
	ax.plot([ϕas_left[end], ϕas_right[end]], [ϕps[end], ϕps[end]], color = "blue", linestyle = "-")


	# Plot binodal lines
	ax.plot(gammas_converter_a(γs, ϕ1s), gammas_converter_p(γs, ϕ1s),
		color = "red", label = L"\mathrm{binodal}")
	ax.plot(gammas_converter_a(γs, ϕ2s), gammas_converter_p(γs, ϕ2s),
		color = "red", label = "_Binodal")

	# Add shading and tie lines based on Pe value

	# find final find gamma
	final_γ = 0.0
	final_ϕ1 = 0.0
	final_ϕ2 = 0.0

	for (γ, ϕ1, ϕ2) in zip(γs, ϕ1s, ϕ2s)
		if (is_stable_value(gamma_converter(γ, ϕ1)...; Pe = Pe)>0) |
		   (is_stable_value(gamma_converter(γ, ϕ2)...; Pe = Pe)>0)
			final_γ = γ
			final_ϕ1 = ϕ1
			final_ϕ2 = ϕ2
			break
		end
	end

	# shading
	tie_line_x = -ϕps*final_γ/(final_γ-1) .+ 1

	ps = collect(0.000001:0.000001:0.4)
	for (γ, ϕ2, ϕ1) in collect(zip(γs, ϕ1s, ϕ2s))[5:10:length(γs)]
		tie_line_x = -ps*γ/(γ-1) .+ 1
		xs = []
		ys = []
		for (x, y) in zip(tie_line_x, ps)
			if (x+y ≤ ϕ1) & (x+y ≥ ϕ2)
				push!(xs, x)
				push!(ys, y)
			end
		end
		if plot_tieline
			ax.plot(xs, ys, color = "grey")
		end
	end

	xs = []
	ys = []
	tie_line_x = -ϕps*final_γ/(final_γ-1) .+ 1
	for (x, y) in zip(tie_line_x, ϕps)
		if (x+y ≤ final_ϕ2) & (x+y ≥ final_ϕ1)
			push!(xs, x)
			push!(ys, y)
		end
	end
	if plot_tieline
		ax.plot([], [], color = "grey", label = L"\mathrm{tie~line}")
	end

	if shading
		ax.fill_betweenx(ϕps, max.(tie_line_x, ϕas_left), ϕas_right, max.(tie_line_x, ϕas_left) .≤ ϕas_right, color = "blue", alpha = 0.3, linewidth = 0)
		ax.fill_betweenx(ϕps, ϕas_left, min.(tie_line_x, ϕas_right), ϕas_left .≤ min.(tie_line_x, ϕas_right), color = "red", alpha = 0.3, linewidth = 0)



		max_ϕa = maximum(ϕas_left)
		ϕas_left, ϕas_right, ϕps, γ_grid, ϕ1_grid, ϕ2_grid = return_spin_from_grid("binodal_1_$(Δϕ)_$(Pe)"; max_ϕa = max_ϕa, Pe = Pe, γ_grid = γs, ϕ1_grid = ϕ1s, ϕ2_grid = ϕ2s, ϕp_grid = gammas_converter_p(γs, ϕ1s))
		ax.fill_betweenx(ϕps, gammas_converter_a(γ_grid, ϕ1_grid), ϕas_left, gammas_converter_a(γ_grid, ϕ1_grid) .≤ ϕas_left, color = "green", alpha = 0.3, linewidth = 0)

		ϕas_left, ϕas_right, ϕps, γ_grid, ϕ1_grid, ϕ2_grid = return_spin_from_grid("binodal_2_$(Δϕ)_$(Pe)"; max_ϕa = max_ϕa, Pe = Pe, γ_grid = γs, ϕ1_grid = ϕ1s, ϕ2_grid = ϕ2s, ϕp_grid = gammas_converter_p(γs, ϕ2s))
		ax.fill_betweenx(ϕps, ϕas_right, gammas_converter_a(γ_grid, ϕ2_grid), gammas_converter_a(γ_grid, ϕ2_grid) .≥ ϕas_right, color = "green", alpha = 0.3, linewidth = 0)

		ax.fill_betweenx(0.0:0.01:0.40, 1.0:(-0.01):0.6, ones(41), color = "grey", alpha = 0.3, linewidth = 0)
	end
	# # plot binodal
	# binod = ax.plot(gammas_converter_a(γs, ϕ1s), gammas_converter_p(γs, ϕ1s), color = "red", label = L"\mathrm{binodal}")
	# ax.plot(gammas_converter_a(γs, ϕ2s), gammas_converter_p(γs, ϕ2s), color = "red", label = "_Bindoal")

	rc("text", usetex = true)
	ax.xaxis.set_tick_params(labelsize = font)
	ax.xaxis.tick_bottom()
	ax.yaxis.set_tick_params(labelsize = font)
	ax.set_xlabel(L"\phi_a", fontsize = font*1.2)
	ax.set_ylabel(L"\phi_p", fontsize = font*1.2)
	ax.tick_params(labelbottom = true, direction = "in")

	# Final axis adjustments
	fig.tight_layout()
	ax.xaxis.set_ticks(0.0:0.2:1.0)
	ax.yaxis.set_ticks(0.0:0.1:0.4)
	ax.axis([0.0, 1.0, 0, 0.4])
	ax.set_aspect((1/(0.4)))
	ax.tick_params(direction = "in", labelsize = font)
	ax.legend(loc = "upper right", fontsize = font, edgecolor = "white")
end

function profile_ρm(axs, font, param, f; lab = "full", col = "black", ls = "-")
	@unpack Δx, Lx, Nx = param

	ax1 = axs[1]
	ax2 = axs[2]
	rho = f[:, 2]+f[:, 1]+f[:, 3];
	rhoa = f[:, 1] + f[:, 2];
	m = f[:, 2] - f[:, 1];
	# Find interface location (minimum of m)
	i = argmin(m);

	# Create x coordinates centered at zero
	xs = (-0.5 .+ (1:Nx)/Nx)

	# Find index closest to x=0
	i_zero = argmin(abs.(xs))

	# Calculate shift needed to move interface to x=0
	shift_amount = i_zero - i

	# Create shifted data
	ρ = circshift(rho, shift_amount);
	m_shifted = circshift(m, shift_amount);


	if lab == "full"
		ax1.plot(xs, ρ;
			color = col, linestyle = ls, label = "\$L = $(Int(Lx))\$")
		ax2.plot(xs, m_shifted;
			color = col, linestyle = ls, label = "\$L = $(Int(Lx))\$")

	else
		ax1.plot(xs, ρ; color = col, linestyle = ls, label = L"L = \infty")

		u   = get_out_u(f, c, Lx)
		Nx  = (length(u)-1)÷2;
		ρ  = u[(0*Nx+1):1:(1*Nx)];
		ρa = u[(1*Nx+1):1:(2*Nx)];

		mouter = zeros(Nx, 1);
		minner = -v0/2*(ds(ρ[1])*ρa[1] + ds(ρ[end])*ρa[end])/2

		ax2.plot(xs, mouter; color = col, linestyle = ls, label = L"L = \infty")
		# ax2.plot(0, minner, "*", color = col)
		ax2.plot(0, 0, ".", color = "white", markersize = 5)
	end


	# ax1.get_xaxis().set_ticks(-0.5:0.25:0.5)# ,fontsize=font)   
	# ax1.get_yaxis().set_ticks(0.4:0.2:1.0)# ,fontsize=font)
	# ax2.get_xaxis().set_ticks(-0.5:0.25:0.5)# ,fontsize=font)
	# ax2.get_yaxis().set_ticks(-0.3:0.1:0.1)# ,fontsize=font)


	# ax1.set_xticklabels([])  # This removes the x-axis tick labels

	# ax2.xaxis.set_tick_params(labelsize=font)
	# ax1.yaxis.set_tick_params(labelsize=font)
	# ax2.yaxis.set_tick_params(labelsize=font)
	# # ax.set_xlabel(L"x", fontsize = 15)
	# #ax.set_ylabel(L"m",fontsize = font, rotation = 90)
	# # ax1.set_aspect((1/(1)))
	# # ax2.set_aspect((1/(1)))
	# ax1.axis([-0.5,0.5,0.4,1.0])
	# # ax2.axis([-0.5,0.5,-0.3,0.1])
	# ax2.axis([-0.5,0.5,-0.3,0.1])
	ax1.tick_params(direction = "in")
	ax1.legend(loc = "upper right", fontsize = font, edgecolor = "white", bbox_to_anchor = (1.02, 1.05))
	ax1.set_ylabel(L"\varrho", fontsize = font)
	ax2.set_xlabel(L"z/L", fontsize = font)
	ax2.set_ylabel(L"\mathrm{m}", fontsize = font)
end

"""
	show_f_reduced(axs, fig, font, param, f; c="?", point=1, Δϕ=0.001, typesol=0)

Shows reduced density profiles and phase space plots.

Arguments:
- axs: Array of two axes for plotting
- fig: Figure object
- font: Font size
- param: Dictionary of parameters
- f: Density profile data
- c: Optional parameter (default "?")
- point: Point index for markers (default 1) 
- Δϕ: Phase space resolution (default 0.001)
- typesol: Type of solution to plot markers for (0-3)
"""
function show_f_reduced(axs, fig, font, param, f; c = "?", point = 1, Δϕ = 0.001, typesol = 0)
	@unpack v0, ϕa, ϕp, Nx = param

	ax1, ax2 = axs[1], axs[2]

	# Plot density profile
	profile_f(ax1, font, param, f)
	ax1.get_xaxis().set_ticks(0:0.25:1)
	ax1.set_xticklabels([L"-0.50", L"-0.25", L"0.0", L"0.25", L"0.50"])
	ax1.set_xlabel(L"z/L", fontsize = font*1.2)

	# Plot phase space
	plot_phase(fig, ax2, param["v0"]; font = font, Δϕ = 0.001)

	ax2.plot(f[:, 1]+f[:, 2], f[:, 3]; color = "black", label = L"(\varrho_a, \varrho_0)")
	ax2.legend(loc = "upper right", fontsize = font, edgecolor = "white")

	# Calculate densities
	ϕp = d3(sum(f)/Nx - sum(f[:, 1:2])/Nx)
	ϕa = d3(sum(f[:, 1:2])/Nx)

	if c > 0
		c = d2(c)
	end

	# Add title with density values
	fig.subplots_adjust(top = 0.85)
	latex_string = latexstring("\$ (\\phi_a, \\phi_p) = ($(ϕa), $(ϕp))\$")
	ax1.set_title(latex_string; fontsize = font)

	# Plot markers based on solution type
	rho = sum(f; dims = 2)[:, 1]
	if typesol == 0
		list_P = []
	elseif typesol == 1
		list_P = zip([Nx÷2, Nx÷2+1], ["P1", "P2"],
			[[0.01, 0.01], [0.01, -0.05]], [[0.01, 0.0], [-0.09, -0.01]])
	elseif typesol == 2
		list_P = zip([Nx÷2, Nx÷2+1, point], ["P2", "P3", "P1"],
			[[0.01, 0.0], [0.01, -0.05], [0.01, -0.02]], [[0.01, 0.0], [-0.08, -0.01], [0.0, 0.005]])
	elseif typesol == 3
		list_P = zip([Nx÷2, Nx÷2+1, point, point+1], ["P3", "P4", "P1", "P2"],
			[[0.01, 0.01], [0.01, -0.05], [0.01, -0.02], [0.00, 0.01]],
			[[0.01, 0.0], [-0.08, 0.0], [0.0, 0.005], [0.01, 0.0]])
	end

	# Add markers and labels
	for (i, label, shift1, shift2) in list_P
		ax1.scatter((i+1)/Nx, rho[i]; c = "none", edgecolors = "black", marker = "o", s = 5, zorder = 3)
		ax1.text((i+1)/Nx+shift1[1], rho[i]+shift1[2], label; fontsize = font, zorder = 3)
		ax2.scatter(f[i, 1]+f[i, 2], f[i, 3]; c = "none", edgecolors = "black", marker = "o", s = 5, zorder = 3)
		ax2.text(f[i, 1]+f[i, 2]+shift2[1], f[i, 3]+shift2[2], label; fontsize = font, zorder = 3)
	end

	ax2.scatter(ϕa, ϕp; c = "black", marker = "^", zorder = 3)
end


function profile_ρ(ax, font, param, f; lab = "full", col = "black", ls = "-")
	@unpack Δx, Lx, Nx = param

	rho = f[:, 2]+f[:, 1]+f[:, 3];
	# Find interface location
	i = argmin(f[:, 2]-f[:, 1]);

	# Create x coordinates centered at zero
	xs = (-0.5 .+ (1:Nx)/Nx)

	# Find index closest to x=0
	i_zero = argmin(abs.(xs))
	println(i_zero)

	# Calculate shift needed to move interface to x=0
	shift_amount = i_zero - i

	println(shift_amount)

	# Create shifted data
	rho_shifted = circshift(rho, shift_amount)

	if lab == "full"
		ax.plot(xs, rho_shifted;
			color = col, linestyle = ls, label = "\$L = $(Int(Lx))\$")
	else
		ax.plot(xs, rho_shifted;
			color = col, linestyle = ls, label = L"L = \infty")
	end


	# normf, c1, dc = f_dot(param, f)
	# latex_string = latexstring("\$ t = $(d2(t)), \\phi_a = $(param["ϕa"]), \\phi_p = $(param["ϕp"]), L_2 = $(d4(normf)), c = $(d4(c1)), {\\dot c} = $(d6(dc))\$")
	# ax.set_title(latex_string, fontsize = font)
	ax.get_xaxis().set_ticks(-0.5:0.25:0.5)# ,fontsize=font)
	# ax.get_yaxis().set_ticks(-0.25:0.25:1.0)#,fontsize=font)
	ax.xaxis.set_tick_params(labelsize = font)
	ax.yaxis.set_tick_params(labelsize = font)
	# ax.set_xlabel(L"x", fontsize = 15)
	#ax.set_ylabel(L"m",fontsize = font, rotation = 90)
	ax.set_aspect((1/(1)))
	ax.axis([-0.5, 0.5, 0.4, 1.0])
	ax.tick_params(direction = "in")
	ax.legend(loc = "upper right", fontsize = font, frameon = false)
	ax.set_xlabel(L"z/L", fontsize = font)
	ax.set_ylabel(L"\rho", fontsize = font)
end


function profile_f(ax, font, param, f)
	@unpack Δx, Lx, Nx = param
	xs = (1/Nx):(1/Nx):1
	ax.plot(xs, f[:, 2]+f[:, 1]+f[:, 3];
		color = "black", linestyle = "-", label = L"\varrho")
	ax.plot(xs, f[:, 1]+f[:, 2];
		color = "red", linestyle = "--", label = L"\varrho_a")
	ax.plot(xs, f[:, 3];
		color = "blue", linestyle = ":", label = L"\varrho_0")
	# ax.plot(xs, (f[:,2]-f[:,1])*Lx; 
	# color = "green", linestyle = "-.", label = L"L_x m")


	# normf, c1, dc = f_dot(param, f)
	# latex_string = latexstring("\$ t = $(d2(t)), \\phi_a = $(param["ϕa"]), \\phi_p = $(param["ϕp"]), L_2 = $(d4(normf)), c = $(d4(c1)), {\\dot c} = $(d6(dc))\$")
	# ax.set_title(latex_string, fontsize = font)
	ax.get_xaxis().set_ticks(0:0.25:1)# ,fontsize=font)
	ax.get_yaxis().set_ticks(-0.25:0.25:1.0)#,fontsize=font)
	ax.xaxis.set_tick_params(labelsize = font)
	ax.yaxis.set_tick_params(labelsize = font)
	# ax.set_xlabel(L"x", fontsize = 15)
	#ax.set_ylabel(L"m",fontsize = font, rotation = 90)
	ax.set_aspect((1/(1)))
	ax.axis([0, 1, 0.0, 1.0])
	ax.tick_params(direction = "in")
	ax.legend(loc = "lower right", fontsize = font, frameon = false)
	ax.set_xlabel(L"{\bar x}", fontsize = font)
end


function plot_phase(fig, ax, Pe; font = 16, Δϕ = 0.01)

	Pes = [Pe]
	axlims = [[0.0, 1.0, 0, 0.4]]
	axs = [ax]
	for (i, (ax, Pe, axlim)) in enumerate(zip(axs, Pes, axlims))

		# load binodal data
		filename = datadir("binodal", "Pe=$(Pe).jld2")
		data = wload(filename)
		@unpack Pe, γs, ϕ1s, ϕ2s = data

		# Set up axis labels and styling
		rc("text", usetex = true)
		ax.xaxis.set_tick_params(labelsize = font)
		ax.xaxis.tick_bottom()
		ax.yaxis.set_tick_params(labelsize = font)
		ax.set_xlabel(L"\phi_a", fontsize = font)
		ax.set_ylabel(L"\phi_p", fontsize = font)
		title = latexstring("\$ \\mathrm{Pe} = $(Pe)\$")
		ax.tick_params(labelbottom = true, direction = "in")

		# Plot spinodal
		ϕas_left, ϕas_right, ϕps, indl, indr = return_spin(; Pe = Pe, Δϕ = Δϕ)
		ax.plot(ϕas_left, ϕps, color = "blue", label = L"\mathrm{spindoal}", linestyle = "-")
		ax.plot(ϕas_right, ϕps, color = "blue", label = "_Spindoal", linestyle = "-")
		ax.plot([ϕas_left[end], ϕas_right[end]], [ϕps[end], ϕps[end]], color = "blue", label = "_Spindoal", linestyle = "-")

		if Pe==7.5
			ax.scatter(ϕas_left[indl], ϕps[indl]; color = "black", marker = "x")
			ax.scatter(ϕas_right[indr], ϕps[indr]; color = "black", marker = "x")
		end

		# phase shading
		if Pe == 5.0
			ax.fill_betweenx(ϕps, ϕas_left, ϕas_right, color = "red", alpha = 0.3, linewidth = 0)
			max_ϕa = maximum(ϕas_left)
			max_ϕp = maximum(ϕps)
			ϕas_left, ϕas_right, ϕps, γ_grid, ϕ1_grid, ϕ2_grid = return_spin_from_grid_real(; max_ϕa = max_ϕa, Pe = Pe, γ_grid = γs, ϕ1_grid = ϕ1s, ϕ2_grid = ϕ2s, ϕp_grid = gammas_converter_p(γs, ϕ1s) .+ 0.00001)
			ax.fill_betweenx(ϕps, gammas_converter_a(γ_grid, ϕ1_grid), ϕas_left, (gammas_converter_a(γ_grid, ϕ1_grid) .≤ ϕas_left), color = "green", alpha = 0.3, linewidth = 0)

			ax.fill_betweenx(ϕps, gammas_converter_a(γ_grid, ϕ1_grid), ϕas_right, (gammas_converter_a(γ_grid, ϕ1_grid) .≥ ϕas_right), color = "green", alpha = 0.3, linewidth = 0)

			ϕas_left, ϕas_right, ϕps, γ_grid, ϕ1_grid, ϕ2_grid = return_spin_from_grid_real(; max_ϕa = max_ϕa, Pe = Pe, γ_grid = γs, ϕ1_grid = ϕ1s, ϕ2_grid = ϕ2s, ϕp_grid = gammas_converter_p(γs, ϕ2s) .+ 0.00001)
			ax.fill_betweenx(ϕps, ϕas_right, gammas_converter_a(γ_grid, ϕ2_grid), gammas_converter_a(γ_grid, ϕ2_grid) .≥ ϕas_right, color = "green", alpha = 0.3, linewidth = 0)

			ax.plot([], [], color = "grey", label = L"\mathrm{tie~line}")

			ps = collect(0.000001:0.000001:0.4)
			for (γ, ϕ2, ϕ1) in collect(zip(γs, ϕ1s, ϕ2s))[2:2:length(γs)]
				tie_line_x = -ps*γ/(γ-1) .+ 1
				xs = []
				ys = []
				for (x, y) in zip(tie_line_x, ps)
					if (x+y ≤ ϕ1)&(x+y ≥ ϕ2)
						push!(xs, x)
						push!(ys, y)
					end
				end
				ax.plot(xs, ys, color = "grey")
			end
		else
			# find final find gamma
			final_γ = 0.0
			final_ϕ1 = 0.0
			final_ϕ2 = 0.0
			for (γ, ϕ1, ϕ2) in zip(γs, ϕ1s, ϕ2s)
				if (is_stable_value(gamma_converter(γ, ϕ1)...; Pe = Pe)>0)|(is_stable_value(gamma_converter(γ, ϕ2)...; Pe = Pe)>0)
					final_γ = γ
					final_ϕ1 = ϕ1
					final_ϕ2 = ϕ2
					break
				end
			end
			# shading
			tie_line_x = -ϕps*final_γ/(final_γ-1) .+ 1
			ax.fill_betweenx(ϕps, max.(tie_line_x, ϕas_left), ϕas_right, max.(tie_line_x, ϕas_left) .≤ ϕas_right, color = "blue", alpha = 0.3, linewidth = 0)
			ax.fill_betweenx(ϕps, ϕas_left, min.(tie_line_x, ϕas_right), ϕas_left .≤ min.(tie_line_x, ϕas_right), color = "red", alpha = 0.3, linewidth = 0)

			ps = collect(0.000001:0.000001:0.4)
			for (γ, ϕ2, ϕ1) in collect(zip(γs, ϕ1s, ϕ2s))[5:10:length(γs)]
				tie_line_x = -ps*γ/(γ-1) .+ 1
				xs = []
				ys = []
				for (x, y) in zip(tie_line_x, ps)
					if (x+y ≤ ϕ1)&(x+y ≥ ϕ2)
						push!(xs, x)
						push!(ys, y)
					end
				end
				ax.plot(xs, ys, color = "grey")
			end

			xs = []
			ys = []
			tie_line_x = -ϕps*final_γ/(final_γ-1) .+ 1
			for (x, y) in zip(tie_line_x, ϕps)
				if (x+y ≤ final_ϕ2)&(x+y ≥ final_ϕ1)
					push!(xs, x)
					push!(ys, y)
				end
			end
			ax.plot([], [], color = "grey", label = L"\mathrm{tie~line}")

			max_ϕa = maximum(ϕas_left)
			ϕas_left, ϕas_right, ϕps, γ_grid, ϕ1_grid, ϕ2_grid = return_spin_from_grid("binodal_1_$(Δϕ)_$(Pe)"; max_ϕa = max_ϕa, Pe = Pe, γ_grid = γs, ϕ1_grid = ϕ1s, ϕ2_grid = ϕ2s, ϕp_grid = gammas_converter_p(γs, ϕ1s))
			ax.fill_betweenx(ϕps, gammas_converter_a(γ_grid, ϕ1_grid), ϕas_left, gammas_converter_a(γ_grid, ϕ1_grid) .≤ ϕas_left, color = "green", alpha = 0.3, linewidth = 0)

			ϕas_left, ϕas_right, ϕps, γ_grid, ϕ1_grid, ϕ2_grid = return_spin_from_grid("binodal_2_$(Δϕ)_$(Pe)"; max_ϕa = max_ϕa, Pe = Pe, γ_grid = γs, ϕ1_grid = ϕ1s, ϕ2_grid = ϕ2s, ϕp_grid = gammas_converter_p(γs, ϕ2s))
			ax.fill_betweenx(ϕps, ϕas_right, gammas_converter_a(γ_grid, ϕ2_grid), gammas_converter_a(γ_grid, ϕ2_grid) .≥ ϕas_right, color = "green", alpha = 0.3, linewidth = 0)

			ax.fill_betweenx(0.0:0.01:0.40, 1.0:(-0.01):0.6, ones(41), color = "grey", alpha = 0.3, linewidth = 0)
		end
		#
		# plot binodal
		binod = ax.plot(gammas_converter_a(γs, ϕ1s), gammas_converter_p(γs, ϕ1s), color = "red", label = L"\mathrm{binodal}")
		ax.plot(gammas_converter_a(γs, ϕ2s), gammas_converter_p(γs, ϕ2s), color = "red", label = "_Bindoal")

		rc("text", usetex = true)
		ax.xaxis.set_tick_params(labelsize = font)
		ax.xaxis.tick_bottom()
		ax.yaxis.set_tick_params(labelsize = font)
		ax.set_xlabel(L"\phi_a", fontsize = font)
		ax.set_ylabel(L"\phi_p", fontsize = font)
		title = latexstring("\$ \\mathrm{Pe} = $(Pe)\$")
		ax.tick_params(labelbottom = true, direction = "in")
		#
	end

	fig.tight_layout()
	ax.xaxis.set_ticks(0.0:0.2:1.0)
	ax.yaxis.set_ticks(0.0:0.1:0.4)
	ax.axis([0.0, 1.0, 0, 0.4])
	ax.tick_params(direction = "in", labelsize = font)
	ax.legend(loc = "upper right", fontsize = font, edgecolor = "white")
end

function show_f(axs, fig, font, param, f; c = "?", point = 1, Δϕ = 0.001)
	@unpack v0, ϕa, ϕp, Nx = param

	ax1 = axs[1]
	ax2 = axs[2]

	profile_f(ax1, font, param, f)

	plot_phase(fig, ax2, param["v0"]; font = font, Δϕ = Δϕ)
	ax2.plot(f[:, 1]+f[:, 2], f[:, 3]; color = "black")

	# check densities
	ϕp = sum(f)/Nx-sum(f[:, 1:2])/Nx
	if ϕp ≈ (param["ϕp"])
		ϕp = param["ϕp"]
	else
		ϕp = d2(ϕp)
	end
	ϕa = sum(f[:, 1:2])/Nx
	if ϕa ≈ (param["ϕa"])
		ϕa = param["ϕa"]
	else
		ϕa = d2(ϕa)
	end
	#

	if c > 0
		c = d2(c)
	end
	normf = sqrt(sum((f[:, 1] .- ϕa/2) .^ 2 + (f[:, 2] .- ϕa/2) .^ 2 + (f[:, 3] .- ϕp) .^ 2)/Nx)

	fig.subplots_adjust(top = 0.85)
	latex_string = latexstring("\$ \\phi_a = $(ϕa), \\phi_p = $(ϕp), \\mathrm{Pe}=$(v0), L_2 = $(d4(normf)), c = $(c)\$")
	fig.suptitle(latex_string; fontsize = font, y = 0.90)

	# plot markers

	rho = sum(f; dims = 2)[:, 1]
	for (i, label, shift1, shift2) in zip([Nx÷2, Nx÷2+1, point], ["B", "C", "A"], [[0.01, 0.01], [0.01, -0.05], [0.00, -0.04]], [[0.01, 0.0], [-0.03, 0.0], [0.0, 0.005]])
		ax1.scatter((i+1)/Nx, rho[i]; c = "none", edgecolors = "black", marker = "o", s = 5, zorder = 3)
		ax1.text((i+1)/Nx+shift1[1], rho[i]+shift1[2], label; fontsize = 9, zorder = 3)
		ax2.scatter(f[i, 1]+f[i, 2], f[i, 3]; c = "none", edgecolors = "black", marker = "o", s = 5, zorder = 3)
		ax2.text(f[i, 1]+f[i, 2]+shift2[1], f[i, 3]+shift2[2], label; fontsize = 9, zorder = 3)
	end

	ax2.scatter(ϕa, ϕp; c = "black", marker = "x", zorder = 3)

	fig_name = "show_f"
	@unpack v0, Lx, Δx = param
	pathname = plotsdir("pm_stretch", "$(fig_name)");
	mkpath(pathname)
	filename = plotsdir("pm_stretch", "$(fig_name)/Lx=$(Lx)_Δx=$(Δx)_Pe=$(v0)_ϕa=$(ϕa)_ϕp=$(ϕp).pdf");
	PyPlot.savefig(filename, dpi = 100, format = "pdf") #bbox_extra_artists=( ldg,)

	display(fig)
end

"""
	create_figure(param)



	Creates a figure showing simulation and PDE results for active particle system. Top row shows kymographs of the simulation and PDE. Bottom row shows the density and magnetization of the simulation and PDE at the final time.
Takes a parameter dictionary as input and returns a figure object.
"""
function create_figure(param)
	# Process data for visualization
	ϵ = 0.1  # Smoothing parameter for local averaging

	images = []
	times  = []
	fts    = []

	# Load simulation and PDE data
	sim_ts, η_saves = load_compress_sim(param)
	pde_ts, f_saves = load_compress_pde(param)

	# Process simulation data
	@unpack N, Nx, N₁, N₂, Lx, Ly = param
	ft_sim = local_average_timeseries(η_saves, ϵ, N, N₁, N₂)
	t_sim_rgb_image = rho_to_rgb(ft_sim)

	# Center simulation pattern
	pk = find_xpeak_ft(sim_ts, ft_sim; time_length = 0.1)
	centre = N₁ ÷ 2 + 1; # place peak at center of domain
	ft_sim = circshift(ft_sim, (0, -pk + centre, 0))
	t_sim_rgb_image = circshift(t_sim_rgb_image, (pk - centre, 0, 0))

	push!(images, t_sim_rgb_image)
	push!(times, sim_ts)
	push!(fts, ft_sim)

	# Process PDE data
	ft_pde = permutedims(reshape(reduce(hcat, f_saves), (Nx, 3, :)), (3, 1, 2))
	t_pde_rgb_image = rho_to_rgb(ft_pde)

	# Center PDE pattern
	pk = find_xpeak_ft(pde_ts, ft_pde; time_length = 0.1)
	centre = Nx ÷ 2 + 1
	ft_pde = circshift(ft_pde, (0, -pk + centre, 0))
	t_pde_rgb_image = circshift(t_pde_rgb_image, (pk - centre, 0, 0))

	push!(images, t_pde_rgb_image)
	push!(times, pde_ts)
	push!(fts, ft_pde)

	# Define plot limits and layout parameters
	rhomin, rhomax = 0.0, 1.0
	mag_lim = 0.8
	font = 12  # Define font size

	# Layout dimensions
	height_1 = 0.08
	width_1 = 0.365
	side_gap_1 = 0.05
	bottom_gap_1 = 0.65
	gap = 0.012

	height_2 = 0.175
	width_2 = height_2
	side_gap_2 = 0.1
	bottom_gap_2 = 0.4
	gap_2 = 0.06

	# Annotation positions
	t_stamp_x = 0.03
	t_stamp_y = 0.05

	# Colorbar positions
	cbar_y_top = bottom_gap_1 + 3*height_1 + 2*gap
	cbar_width = 0.1
	cbar_y_bot = bottom_gap_1 + gap + height_1
	cbar_x = 0.885

	# Initialize figure
	rc("text", usetex = true)
	fig = plt.figure(figsize = (10, 10))

	# Create base axes for the figure
	ax = fig.add_axes([0, 0, 1, 1], visible = true)
	ax.xaxis.set_ticks([])
	ax.yaxis.set_ticks([])
	for spine in ["top", "right", "bottom", "left"]
		ax.spines[spine].set_visible(false)
	end

	@unpack T = param

	# Add image plots
	for (i, (rgb_image, ts)) in enumerate(zip(images, times))
		ax = fig.add_axes([side_gap_1+(i-1)*(side_gap_1+width_1), bottom_gap_1, width_1, height_1])
		t_end = ts[end]
		t_start = ts[1]
		ax.imshow(rgb_image; extent = [t_start, t_end, 0, Lx], interpolation = "bilinear")
		ax.get_yaxis().set_ticks(0:1.0:Lx)
		ax.set_yticklabels(["0", "1", "2"])
		# ax.get_xaxis().set_ticks([])
		ax.axis([0, T, 0, 2])
		ax.set_aspect((T/2)*(height_1/width_1))
		ax.set_ylabel(L"x", fontsize = font, rotation = 90)
		# ax.get_xaxis().set_ticks(0:round(0.25*round(T); digits = 0):round(T;digits = 0 ))
		ax.set_xlabel(L"t", fontsize = font)
		ax.tick_params(labelbottom = true, direction = "in")
	end

	# Add matrix plots
	for (i, (ts, ft)) in enumerate(zip(times, fts))
		global im1, im2
		# Density plot
		ax = fig.add_axes([side_gap_1+(i-1)*(side_gap_1+width_1), bottom_gap_1+2*height_1+2*gap, width_1, height_1])
		t_end = ts[end]
		t_start = ts[1]
		_, Nx, _ = size(ft)

		colmap = PyPlot.plt.cm.viridis
		norm1 = matplotlib.colors.Normalize(vmin = rhomin, vmax = rhomax)
		im1 = ax.matshow((ft[:, Nx:-1:1, 1]+ft[:, Nx:-1:1, 2]+ft[:, Nx:-1:1, 3])';
			norm = norm1, cmap = colmap, extent = [t_start, t_end, 0, Lx])

		ax.xaxis.set_ticks([])
		ax.xaxis.tick_bottom()
		ax.get_yaxis().set_ticks(0:1.0:Lx)
		ax.set_yticklabels(["0", "1", "2"])
		ax.axis([0, T, 0, 2])
		ax.set_aspect((T/2)*(height_1/width_1))
		ax.set_ylabel(L"x", fontsize = font, rotation = 90)
		# ax.get_xaxis().set_ticks(0:round(0.25*round(T); digits = 0):round(T;digits = 0 ))
		ax.tick_params(labelbottom = false, direction = "in")

		# Magnetization plot
		ax = fig.add_axes([side_gap_1+(i-1)*(side_gap_1+width_1), bottom_gap_1+height_1+gap, width_1, height_1])
		colmap = PyPlot.plt.cm.PRGn
		norm1 = matplotlib.colors.Normalize(vmin = -mag_lim, vmax = mag_lim)
		im2 = ax.matshow((ft[:, Nx:-1:1, 2]-ft[:, Nx:-1:1, 1])';
			norm = norm1, cmap = colmap, extent = [t_start, t_end, 0, Lx])

		# ax.get_xaxis().set_ticks([])
		ax.get_yaxis().set_ticks(0:1.0:Lx)
		ax.set_yticklabels(["0", "", "2"])
		ax.xaxis.tick_bottom()
		ax.axis([0, T, 0, 2])
		ax.set_aspect((T/2)*(height_1/width_1))
		ax.set_ylabel(L"x", fontsize = font, rotation = 90)
		# ax.get_xaxis().set_ticks(0:round(0.25*round(T); digits = 0):round(T;digits = 0 ))
		ax.tick_params(labelbottom = false, direction = "in")
	end

	# Add colorbars
	cbar_ax = fig.add_axes([cbar_x, bottom_gap_1, height_1, height_1])
	Δ = 0.001
	cbar_f = [x*(x+y≤1)*(i!=3)/2 + y*(x+y≤1)*(i==3) for x in Δ:Δ:1, y in Δ:Δ:1, i in 1:3]
	rgb_image = rho_to_rgb(cbar_f)

	ax = cbar_ax
	ax.imshow(rgb_image; extent = [0, 1, 0, 1])
	ax.spines["top"].set_visible(false)
	ax.spines["right"].set_visible(false)
	ax.get_xaxis().set_ticks(0:0.5:1)
	ax.get_yaxis().set_ticks(0:0.5:1)
	ax.set_xlabel(L"\rho_a", fontsize = font)
	ax.set_ylabel(L"\rho_0", fontsize = font, rotation = 90)
	ax.tick_params(direction = "in")

	# Add density colorbar
	rho_cbar_ax = fig.add_axes([cbar_x, cbar_y_bot+gap+height_1, 0.025, height_1])
	rho_cbar = fig.colorbar(im1, cax = rho_cbar_ax)
	rho_cbar_ax.set_ylabel(L"\rho", fontsize = font, rotation = 90)
	rho_cbar.set_ticks(rhomin:0.5:rhomax)
	rho_cbar_ax.yaxis.set_ticks_position("right")
	rho_cbar_ax.tick_params(direction = "in")

	# Add magnetization colorbar
	mag_cbar_ax = fig.add_axes([cbar_x, cbar_y_bot, 0.025, height_1])
	mag_cbar = fig.colorbar(im2, cax = mag_cbar_ax)
	mag_cbar.set_ticks((-mag_lim):0.8:mag_lim)
	mag_cbar_ax.tick_params(direction = "in")
	mag_cbar_ax.yaxis.set_ticks_position("right")
	mag_cbar_ax.set_ylabel(L"m", fontsize = font, rotation = 90)

	# Add final time plots
	# Total density (ρ)
	ax = fig.add_axes([side_gap_2, bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 2] + fts[1][end, :, 1] + fts[1][end, :, 3]
	pde_mag = fts[2][end, :, 2] + fts[2][end, :, 1] + fts[2][end, :, 3]

	@unpack N, Δx = param
	ax.plot(Δx:Δx:Lx, pde_mag; color = "black", label = "PDE")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red", label = "Simulation")

	ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks(rhomin:0.5:rhomax)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"\rho", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(rhomax-rhomin)))
	ax.axis([0, Lx, rhomin, rhomax])
	ax.tick_params(direction = "in")

	# Active density (ρₐ)
	ax = fig.add_axes([side_gap_2+gap_2+width_2, bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 2] + fts[1][end, :, 1]
	pde_mag = fts[2][end, :, 2] + fts[2][end, :, 1]

	ax.plot(Δx:Δx:Lx, pde_mag; color = "black")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red")

	ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks(rhomin:0.5:rhomax)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"\rho_a", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(rhomax-rhomin)))
	ax.axis([0, Lx, rhomin, rhomax])
	ax.tick_params(direction = "in")

	# Passive density (ρ₀)
	ax = fig.add_axes([side_gap_2+2*(gap_2+width_2), bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 3]
	pde_mag = fts[2][end, :, 3]

	ax.plot(Δx:Δx:Lx, pde_mag; color = "black")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red")

	ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks(rhomin:0.5:rhomax)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"\rho_0", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(rhomax-rhomin)))
	ax.axis([0, Lx, rhomin, rhomax])
	ax.tick_params(direction = "in")

	# Magnetization (m)
	ax = fig.add_axes([side_gap_2+3*(gap_2+width_2), bottom_gap_2, width_2, height_2])

	sim_mag = fts[1][end, :, 2] - fts[1][end, :, 1]
	pde_mag = fts[2][end, :, 2] - fts[2][end, :, 1]

	ax.plot(Δx:Δx:Lx, pde_mag; color = "black")
	ax.plot((1/N):(1/N):Lx, sim_mag; color = "red")

	ax.get_xaxis().set_ticks(0:0.5:Lx)
	ax.get_yaxis().set_ticks((-mag_lim):0.4:mag_lim)
	ax.set_xlabel(L"x", fontsize = font)
	ax.set_ylabel(L"m", fontsize = font, rotation = 90)
	ax.set_aspect((Lx/(2*mag_lim)))
	ax.axis([0, Lx, -mag_lim, mag_lim])
	ax.tick_params(direction = "in")

	# Add time stamp and labels
	latex_annotation = latexstring("\$ t = $(round(times[1][end];digits = 1))\$")
	ax.annotate(latex_annotation, (t_stamp_x, bottom_gap_2+t_stamp_y),
		xycoords = "figure fraction", rotation = 90, fontsize = font)

	# Add subplot labels
	ax.annotate(L"(a)", (side_gap_1 - 0.03, bottom_gap_1+3*height_1+2*gap),
		xycoords = "figure fraction", rotation = 0, fontsize = font)
	ax.annotate(L"(b)", (2*side_gap_1+width_1 - 0.03, bottom_gap_1+3*height_1+2*gap),
		xycoords = "figure fraction", rotation = 0, fontsize = font)
	ax.annotate(L"(c)", (side_gap_1, bottom_gap_2+height_2+0.01),
		xycoords = "figure fraction", rotation = 0, fontsize = font)

	return fig;
end


function rho_to_rgb(f)
	Nx, Ny, k = size(f)
	rgb_image = ones(Ny, Nx, 3)

	rgb_image[:, :, 3] = -(f[:, Ny:-1:1, 1]' + f[:, Ny:-1:1, 2]') .^ 2 .+ 1
	rgb_image[:, :, 1] = -(f[:, Ny:-1:1, 3]') .+ 1
	rgb_image[:, :, 2] = -(f[:, Ny:-1:1, 1]' + f[:, Ny:-1:1, 2]' - f[:, Ny:-1:1, 3]') .^ 2 .+ 1

	return rgb_image
end

function find_xpeak_ft(ts, ft; time_length = 0.1)
	index_end = length(ts)
	indebottom_gap_2 = index_end - length([t for t in ts if t>ts[end]-time_length])

	av_rho = sum(ft[indebottom_gap_2:1:index_end, :, :]; dims = (1, 3))[1, :, 1]
	N = length(av_rho)
	M = 3*N
	av_rho = reduce(vcat, [av_rho, av_rho, av_rho])

	smooth_rho = KernelDensitySJ.smooth(collect(1.0:M), av_rho, Int64(N/10 ÷ 1), collect(1.0:M))
	pks, vals = findmaxima(smooth_rho)
	pks = [x for x in pks if (x>M/6)&(x<5*M/6)]
	pks, proms = peakproms(pks, smooth_rho)
	pk = pks[argmax(proms)]
	return pk
end

# 
println("v5.0")
