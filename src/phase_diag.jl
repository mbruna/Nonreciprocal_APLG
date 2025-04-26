# Functions to perform phase diagram calculations of the active particle lattice gas

using QuadGK, ForwardDiff, Roots, LinearAlgebra

# Define functions 

rtol = 1e-14
tol = 1e-14
d3(x) = round(x; digits = 3);


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

ds(x) = self_diff(x);
dsp(x) = self_diff_prime(x);

F(x) = (1-x)*(ds(x) + dsp(x)*x)
pspin(x; Pe = 10.0) = (-2-F(x)*Pe^2)/(ds(x)-dsp(x)*(1-x))/Pe^2

R_prime(ρ) = 1/ds(ρ)/(1-ρ)^3


########################
function R(ρ::Float64; rtol = rtol, tol = tol)
	α::Float64 = π/2 - 1;
	c::Float64 = -1/(α*(2*α-1)/(2*α+1) - α + 1);
	f(x) = R_prime(x) - c/(1-x)^4;
	return quadgk(f, 0.0, ρ; rtol = rtol, atol = tol)[1] + c/3/(1-ρ)^3 - c/3
end

function R(x::ForwardDiff.Dual{T, V, N}; rtol = rtol, tol = tol) where {T, V, N}
	return ForwardDiff.Dual{T}(R(x.value; rtol = rtol, tol = tol), R_prime(x.value)*ForwardDiff.partials(x))
end


κ(ρ, Pe) = self_diff(ρ) / Pe / (1 - ρ)
Λ(ρ, Pe) = -2 * self_diff(ρ) / Pe / (1 - ρ)^2
g0(ρ; Pe = 10, γ = 1, atol = 1) = Pe*(1-γ*(1 - ρ))*self_diff(ρ) - 2*log(1 - ρ)/Pe
dg0(ρ; Pe = 10, γ = 1) = Pe*(γ)*self_diff(ρ) + Pe*(1-γ*(1 - ρ))*self_diff_prime(ρ) + 2/(1 - ρ)/Pe


g0(x::ForwardDiff.Dual{T, V, N}; Pe = 10, γ::ForwardDiff.Dual{T, V, N} = 1, atol = 1) where {T, V, N} = ForwardDiff.Dual{T}(g0(x.value; Pe = Pe, γ = γ), dg0(x.value; Pe = Pe, γ = γ)*ForwardDiff.partials(x))


dΦ_dρ(ρ, Pe, γ) = g0(ρ; Pe = Pe, γ = γ)*R_prime(ρ)

function ΦoR(x::Float64; Pe = 10, γ = 3.4, rtol = rtol, tol = tol) # Φ as a funciton of ρ ie composed with R
	α::Float64 = π/2 - 1;
	c::Float64 = 2/Pe/(α*(2*α-1)/(2*α+1) - α + 1);
	f(x) = dΦ_dρ(x, Pe, γ) - c*log(1 - x)/(1-x)^4;
	return quadgk(f, 0.0, x, rtol = rtol, atol = tol)[1]+c*(1+3*log(1-x))/9/(1-x)^3-c*(1+3*log(1-0.95))/9/(1-0.95)^3
end

function ΦoR(x::ForwardDiff.Dual{T, V, N}; Pe = 10, γ = 3.4, rtol = rtol, tol = tol) where {T, V, N} # Φ as a funciton of ρ ie composed with R
	return return ForwardDiff.Dual{T}(ΦoR(x.value), dΦ_dρ(x.value, Pe, γ)*ForwardDiff.partials(x))
end

function h0(x::ForwardDiff.Dual{T, V, N}; Pe = 10, γ = 1, atol = 1e-12) where {T, V, N}
	return ForwardDiff.Dual{T}(g0(x.value; Pe = Pe, γ = γ)*R.(x.value; tol = atol)-ΦoR.(x.value; Pe = Pe, γ = γ, tol = atol), (dg0(x.value; Pe = Pe, γ = γ)*R.(x.value; tol = atol))*ForwardDiff.partials(x))
end

function h0(ρ::Float64; Pe = 10, γ = 1, atol = 1e-12)
	return g0(ρ; Pe = Pe, γ = γ)*R.(ρ; tol = atol)-ΦoR.(ρ; Pe = Pe, γ = γ, tol = atol)
end

#### 
# Find solution

function find_G_density_no_error(ϕ1; Pe = 5, γ = 1, limits = (0.8, 0.9999), fail = 0.0, atol = 1e-12)
	try
		f(x) = g0(x; Pe = Pe, γ = γ, atol = atol)-g0(ϕ1; Pe = Pe, γ = γ, atol = atol)
		ϕ2 = find_zero(f, limits)
		return ϕ2
	catch
		return fail
	end
end

function find_H_density_no_error(ϕ1; Pe = 5, γ = 1, limits = (0.8, 0.9999), fail = 0.0, atol = 1e-12)
	try
		f(x) = h0(x; Pe = Pe, γ = γ, atol = atol)-h0(ϕ1; Pe = Pe, γ = γ, atol = atol)
		ϕ2 = find_zero(f, limits)
		return ϕ2
	catch
		return fail
	end
end

#

function g_tunring_points(; Pe = 5, γ = 3.52, initial_Δ = 1e-4, atol = atol)
	x = 0.0:(initial_Δ):(1-initial_Δ)
	gmin, gmax = (0.0, 1.0)
	DG = dg0.(x; Pe = Pe, γ = γ)
	DG_min = minimum(DG)
	gmin = 1.0
	gmax = 0.0
	if DG_min > 0
		return println("no turn in g0")
	else
		i = argmin(DG)
		x0 = x[i]
		f(x) = dg0(x; Pe = Pe, γ = γ)
		gmax = find_zero(f, (0.0, x0))
		gmin = find_zero(f, (x0, 1.0))
	end
	return gmin, gmax
end

function h_maximum(; Pe = 5.0, γ = 1.0, initial_Δ = 1e-4, atol = atol)
	x = 0.0:(initial_Δ):(1-initial_Δ)
	return x[argmax(h0.(x; Pe = Pe, γ = γ, atol = atol))]
end

function h_minimum(hmax; Pe = 5.0, γ = 1.0, initial_Δ = 1e-4, atol = atol)
	X = 0.0:(initial_Δ):(1-initial_Δ)
	x = [x for x in X if x < hmax]
	return x[argmin(h0.(x; Pe = Pe, γ = γ, atol = atol))]
end

function initial_intervals(; Pe = 5, γ = 1, rho_max = (1-10e-25), initial_Δ = 1e-4, atol = atol)
	local lower_limits, upper_limits
	gmin, gmax = g_tunring_points(; Pe = Pe, γ = γ, initial_Δ = initial_Δ, atol = atol)
	lower_1 = find_G_density_no_error(gmin; Pe = Pe, γ = γ, limits = (0.0, gmax), fail = 0.0, atol = atol)
	upper_2 = find_G_density_no_error(gmax; Pe = Pe, γ = γ, limits = (gmin, rho_max), fail = rho_max, atol = atol)

	#h_max = h_minimum(;Pe = Pe, γ = γ, initial_Δ = initial_Δ, atol = atol)

	lower_limits = (lower_1, gmax)
	upper_limits = (gmin, upper_2)
	return lower_limits, upper_limits
end

function restrict_intervals_h(lower_limits, upper_limits; Pe = 5, γ = 1, atol = atol)
	upper_1 = find_H_density_no_error(lower_limits[1]; Pe = Pe, γ = γ, limits = upper_limits, fail = upper_limits[1], atol = atol)
	upper_2 = find_H_density_no_error(lower_limits[2]; Pe = Pe, γ = γ, limits = upper_limits, fail = upper_limits[2], atol = atol)
	upper_limits = (min(upper_1, upper_2), max(upper_1, upper_2))
	return lower_limits, upper_limits
end

function restrict_intervals_g(lower_limits, upper_limits; Pe = 5, γ = 1, atol = atol)
	lower_1 = find_G_density_no_error(upper_limits[1]; Pe = Pe, γ = γ, limits = lower_limits, fail = lower_limits[1], atol = atol)
	lower_2 = find_G_density_no_error(upper_limits[2]; Pe = Pe, γ = γ, limits = lower_limits, fail = lower_limits[2], atol = atol)
	lower_limits = (min(lower_1, lower_2), max(lower_1, lower_2))
	return lower_limits, upper_limits
end

function interval_size(lower_limits, upper_limits)
	return abs(lower_limits[2]-lower_limits[1])+abs(upper_limits[2]-upper_limits[1])
end

function colapse_sol_interval(; Pe = 5, γ = 1.0, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 20, tol = 1e-8, atol = 1e-12)
	lower_limits, upper_limits = initial_intervals(; Pe = Pe, γ = γ, rho_max = rho_max, initial_Δ = initial_Δ, atol = atol)
	precision = interval_size(lower_limits, upper_limits)
	i = 0
	while (i<max_iter)&(precision>tol)
		lower_limits, upper_limits = restrict_intervals_h(lower_limits, upper_limits; Pe = Pe, γ = γ, atol = atol)
		lower_limits, upper_limits = restrict_intervals_g(lower_limits, upper_limits; Pe = Pe, γ = γ, atol = atol)
		precision = interval_size(lower_limits, upper_limits)
		i+=1
	end
	if precision ≤ tol
		find_sol = true
	else
		find_sol = false
	end

	return find_sol, lower_limits, upper_limits
end

#

function dg0_minimum(; Pe = 5.0, γ = 1.0, initial_Δ = 1e-4)
	x = 0.0:(initial_Δ):(1-initial_Δ)
	DG = dg0.(x; Pe = Pe, γ = γ)
	return minimum(DG)
end

function dg0_argmin(; Pe = 5.0, γ = 1.0, initial_Δ = 1e-4)
	x = 0.0:(initial_Δ):(1-initial_Δ)
	DG = dg0.(x; Pe = Pe, γ = γ)
	return argmin(DG)*initial_Δ
end

function find_gamma_limit(; Pe = 5.0, initial_Δ = 1e-5, γ_max = 100.0)
	f(x) = dg0_minimum(; Pe = Pe, γ = x, initial_Δ = initial_Δ)
	return find_zero(f, (1.0, γ_max))
end

function find_rho_limit(; Pe = 5.0, initial_Δ = 1e-5, γ_max = 100.0)
	γ = find_gamma_limit(; Pe = Pe, initial_Δ = initial_Δ, γ_max = γ_max)
	return dg0_argmin(; Pe = Pe, γ = γ, initial_Δ = initial_Δ)
end

###

function chi_converter(γ, ϕ)
	ϕa = 1 - γ*(1 - ϕ)
	return ϕa/ϕ
end

function chis_converter(γs, ϕs)
	map(γs, ϕs) do γ, ϕ
		chi_converter(γ, ϕ)
	end
end

function gamma_converter(γ, ϕ)
	ϕa = 1 - γ*(1 - ϕ)
	ϕp = ϕ-ϕa
	return ϕa, ϕp
end

function gammas_converter_a(γs, ϕs)
	map(γs, ϕs) do γ, ϕ
		gamma_converter(γ, ϕ)[1]
	end
end

function gammas_converter_p(γs, ϕs)
	map(γs, ϕs) do γ, ϕ
		gamma_converter(γ, ϕ)[2]
	end
end

###

function α(ϕa, ϕp)
	ϕ  = ϕa + ϕp
	ϕ0 = 1 - ϕ
	ds  = self_diff(ϕ)
	dsp = self_diff_prime(ϕ)
	return ϕa*ϕ0*(ds+dsp*ϕ)+ds^2*ϕp
end

function β(ϕa, ϕp)
	ϕ  = ϕa + ϕp
	ϕ0 = 1 - ϕ
	ds  = self_diff(ϕ)
	dsp = self_diff_prime(ϕ)
	return 2*(1-ds)*ϕ
end

function is_complex_value(ϕa, ϕp; Pe = 10.0)
	ϕ = ϕa + ϕp
	ϕ0 = 1 - ϕ
	ds = self_diff(ϕ)
	dsp = self_diff_prime(ϕ)
	expr2 = (β(ϕa, ϕp) + α(ϕa, ϕp)*Pe^2)
	return -(expr2^2 - 4*ϕp*β(ϕa, ϕp)*ds^2*Pe^2)
end

function u1dw(p, pa, Pe, w, ds, dds)
	return (
		-6*(2*p^2*w + 4*ds*p^2*w)*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^2 +
		9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
		(4*p^4*w + 4*ds*p^4*w + 2*ds^2*p^4*Pe^2*w + 2*ds*p^3*pa*Pe^2*w - 2*ds^2*p^3*pa*Pe^2*w +
		 2*dds*p^4*pa*Pe^2*w - 2*ds*p^4*pa*Pe^2*w - 2*dds*p^5*pa*Pe^2*w + 8*ds*p^4*w^3 + 4*ds^2*p^4*w^3) +
		9*(2*p^2*w + 4*ds*p^2*w)*(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 -
								  ds^2*p^3*pa*Pe^2*w^2 + dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 +
								  2*ds*p^4*w^4 + ds^2*p^4*w^4) - 27*(8*ds*p^6*w^3 + 4*ds^2*p^6*Pe^2*w^3 + 4*dds*ds*p^6*pa*Pe^2*w^3 -
			4*ds^2*p^6*pa*Pe^2*w^3 - 4*dds*ds*p^7*pa*Pe^2*w^3 + 6*ds^2*p^6*w^5)
	)/54 +
		   (
		(p^9 * (2*(-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2)*w -
				4*(-1 + ds)^2*p*w^3) *
		 (-4*p + (-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2)*w^2 -
		  (-1 + ds)^2*p*w^4)^2)/243 +
		(
			(
				6*(2*p^2*w + 4*ds*p^2*w)*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^2 -
				9*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)*(4*p^4*w + 4*ds*p^4*w + 2*ds^2*p^4*Pe^2*w +
													2*ds*p^3*pa*Pe^2*w - 2*ds^2*p^3*pa*Pe^2*w + 2*dds*p^4*pa*Pe^2*w - 2*ds*p^4*pa*Pe^2*w -
													2*dds*p^5*pa*Pe^2*w + 8*ds*p^4*w^3 + 4*ds^2*p^4*w^3) -
				9 * (2*p^2*w + 4*ds*p^2*w) *
				(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
				 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) +
				27*(8*ds*p^6*w^3 + 4*ds^2*p^6*Pe^2*w^3 + 4*dds*ds*p^6*pa*Pe^2*w^3 - 4*ds^2*p^6*pa*Pe^2*w^3 -
					4*dds*ds*p^7*pa*Pe^2*w^3 + 6*ds^2*p^6*w^5)
			)*(
				2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 -
				9*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)*(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 +
													ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 + dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 -
													dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) +
				27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 - ds^2*p^6*pa*Pe^2*w^4 -
					dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
			)
		)/1458
	) /
		   (
		2*sqrt(
			(p^9*(-4*p + (-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2) *
						 w^2 - (-1 + ds)^2*p*w^4)^3)/729 +
			(
				2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 -
				9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
				(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
				 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 +
				 ds^2*p^4*w^4) + 27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 -
					ds^2*p^6*pa*Pe^2*w^4 - dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
			)^2/2916,
		)
	);
end

function u2dw(p, pa, Pe, w, ds, dds)
	return (
		-6*(2*p^2*w + 4*ds*p^2*w)*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^2 +
		9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
		(4*p^4*w + 4*ds*p^4*w + 2*ds^2*p^4*Pe^2*w + 2*ds*p^3*pa*Pe^2*w - 2*ds^2*p^3*pa*Pe^2*w +
		 2*dds*p^4*pa*Pe^2*w - 2*ds*p^4*pa*Pe^2*w - 2*dds*p^5*pa*Pe^2*w + 8*ds*p^4*w^3 + 4*ds^2*p^4*w^3) +
		9*(2*p^2*w + 4*ds*p^2*w)*(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 -
								  ds^2*p^3*pa*Pe^2*w^2 + dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 +
								  2*ds*p^4*w^4 + ds^2*p^4*w^4) - 27*(8*ds*p^6*w^3 + 4*ds^2*p^6*Pe^2*w^3 + 4*dds*ds*p^6*pa*Pe^2*w^3 -
			4*ds^2*p^6*pa*Pe^2*w^3 - 4*dds*ds*p^7*pa*Pe^2*w^3 + 6*ds^2*p^6*w^5)
	)/54 -
		   (
		(p^9 * (2*(-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2)*w -
				4*(-1 + ds)^2*p*w^3) *
		 (-4*p + (-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2)*w^2 -
		  (-1 + ds)^2*p*w^4)^2)/243 +
		(
			(
				6*(2*p^2*w + 4*ds*p^2*w)*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^2 -
				9*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)*(4*p^4*w + 4*ds*p^4*w + 2*ds^2*p^4*Pe^2*w +
													2*ds*p^3*pa*Pe^2*w - 2*ds^2*p^3*pa*Pe^2*w + 2*dds*p^4*pa*Pe^2*w - 2*ds*p^4*pa*Pe^2*w -
													2*dds*p^5*pa*Pe^2*w + 8*ds*p^4*w^3 + 4*ds^2*p^4*w^3) -
				9 * (2*p^2*w + 4*ds*p^2*w) *
				(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
				 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) +
				27*(8*ds*p^6*w^3 + 4*ds^2*p^6*Pe^2*w^3 + 4*dds*ds*p^6*pa*Pe^2*w^3 - 4*ds^2*p^6*pa*Pe^2*w^3 -
					4*dds*ds*p^7*pa*Pe^2*w^3 + 6*ds^2*p^6*w^5)
			)*(
				2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 -
				9*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)*(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 +
													ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 + dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 -
													dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) +
				27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 - ds^2*p^6*pa*Pe^2*w^4 -
					dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
			)
		)/1458
	) /
		   (
		2*sqrt(
			(p^9*(-4*p + (-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2) *
						 w^2 - (-1 + ds)^2*p*w^4)^3)/729 +
			(
				2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 -
				9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
				(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
				 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 +
				 ds^2*p^4*w^4) + 27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 -
					ds^2*p^6*pa*Pe^2*w^4 - dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
			)^2/2916,
		)
	);
end

function shiftdw(p, pa, Pe, w, ds, dds)
	return (-2*p^2*w - 4*ds*p^2*w)/3
end

function u1(p, pa, Pe, w, ds, dds)
	return (
		-2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 +
		9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
		(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
		 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) -
		27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 - ds^2*p^6*pa*Pe^2*w^4 -
			dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
	)/54 +
		   sqrt(
		(p^9*(-4*p + (-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2)*w^2 -
			  (-1 + ds)^2*p*w^4)^3)/729 +
		(
			2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 -
			9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
			(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
			 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) +
			27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 - ds^2*p^6*pa*Pe^2*w^4 -
				dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
		)^2/2916,
	)
end

function u2(p, pa, Pe, w, ds, dds)
	return (
		-2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 +
		9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
		(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
		 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) -
		27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 - ds^2*p^6*pa*Pe^2*w^4 -
			dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
	)/54 -
		   sqrt(
		(p^9*(-4*p + (-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2)*w^2 -
			  (-1 + ds)^2*p*w^4)^3)/729 +
		(
			2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 -
			9 * (2*p^2 + p^2*w^2 + 2*ds*p^2*w^2) *
			(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 + ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 +
			 dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 - dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4) +
			27*(2*ds*p^6*w^4 + ds^2*p^6*Pe^2*w^4 + dds*ds*p^6*pa*Pe^2*w^4 - ds^2*p^6*pa*Pe^2*w^4 -
				dds*ds*p^7*pa*Pe^2*w^4 + ds^2*p^6*w^6)
		)^2/2916,
	);
end


function shift_fn(p, pa, Pe, w, ds, dds)
	return (-2*p^2 - p^2*w^2 - 2*ds*p^2*w^2)/3
end

function re_eigen_f(ϕ, ϕa, v0, w, ds, dsp)
	return (-(cbrt(u1(ϕ, ϕa, v0, w, ds, dsp)) + cbrt(u2(ϕ, ϕa, v0, w, ds, dsp)))/2 +
			shift_fn(ϕ, ϕa, v0, w, ds, dsp))/ϕ^2
end

function discriminant_f(p, pa, Pe, w, ds, dds)
	return (p^9*(-4*p + (-2*(-1 + ds)*p + 3*(ds^2*p - (ds^2 + ds*(-1 + p) + dds*(-1 + p)*p)*pa)*Pe^2)*w^2 -
				 (-1 + ds)^2*p*w^4)^3)/729 +
		   (
		2*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)^3 -
		9*(2*p^2 + p^2*w^2 + 2*ds*p^2*w^2)*(2*p^4*w^2 + 2*ds*p^4*w^2 + ds^2*p^4*Pe^2*w^2 +
											ds*p^3*pa*Pe^2*w^2 - ds^2*p^3*pa*Pe^2*w^2 + dds*p^4*pa*Pe^2*w^2 - ds*p^4*pa*Pe^2*w^2 -
											dds*p^5*pa*Pe^2*w^2 + 2*ds*p^4*w^4 + ds^2*p^4*w^4)
	)^2/2916
end

function Dw(ϕ, ϕa, v0, w, ds, dsp)
	# if discriminant_f(ϕ, ϕa, v0, w, ds, dsp) ≥ 0
	#     return (-(  u1dw(ϕ, ϕa, v0, w, ds, dsp)*cbrt(u1(ϕ, ϕa, v0, w, ds, dsp))^(-2)+  u2dw(ϕ, ϕa, v0, w, ds, dsp)*cbrt(u2(ϕ, ϕa, v0, w, ds, dsp))^(-2)    )/6 + shiftdw(ϕ, ϕa, v0, w, ds, dsp))/ϕ^2
	# else
	vars = complex.([ϕ, ϕa, v0, w, ds, dsp])
	return real(((u1dw(vars...)*(u1(vars...))^(-2/3)+u2dw(vars...)*(u2(vars...))^(-2/3))/3 + shiftdw(vars...))/ϕ^2)
	# end
end

function Dw2(ϕ, ϕa, v0, w, ds, dsp)
	# if discriminant_f(ϕ, ϕa, v0, w, ds, dsp) ≥ 0
	#     return ((  u1dw(ϕ, ϕa, v0, w, ds, dsp)*cbrt(u1(ϕ, ϕa, v0, w, ds, dsp))^(-2)+  u2dw(ϕ, ϕa, v0, w, ds, dsp)*cbrt(u2(ϕ, ϕa, v0, w, ds, dsp))^(-2)    )/3 + shiftdw(ϕ, ϕa, v0, w, ds, dsp))/ϕ^2
	# else
	vars = complex.([ϕ, ϕa, v0, w, ds, dsp])
	z = (-1+im*sqrt(3))/2
	return real(((z*u1dw(vars...)*(u1(vars...))^(-2/3)+z^2*u2dw(vars...)*(u2(vars...))^(-2/3))/3 + shiftdw(vars...))/ϕ^2)
	# end
end

function Dw3(ϕ, ϕa, v0, w, ds, dsp)
	# if discriminant_f(ϕ, ϕa, v0, w, ds, dsp) ≥ 0
	#     return ((  u1dw(ϕ, ϕa, v0, w, ds, dsp)*cbrt(u1(ϕ, ϕa, v0, w, ds, dsp))^(-2)+  u2dw(ϕ, ϕa, v0, w, ds, dsp)*cbrt(u2(ϕ, ϕa, v0, w, ds, dsp))^(-2)    )/3 + shiftdw(ϕ, ϕa, v0, w, ds, dsp))/ϕ^2
	# elsed
	vars = complex.([ϕ, ϕa, v0, w, ds, dsp])
	z = (-1+im*sqrt(3))/2
	return real(((z^2*u1dw(vars...)*(u1(vars...))^(-2/3)+z*u2dw(vars...)*(u2(vars...))^(-2/3))/3 + shiftdw(vars...))/ϕ^2)
	# end
end

function relambw(w; ϕa = 0.1, ϕ = 0.7, v = 7.5)
	ω = w;
	Pe = v;
	ϕ0 = 1 - ϕ;
	DS = self_diff(ϕ);
	DSP = self_diff_prime(ϕ);

	DD = (1-DS)/ϕ
	s = DD - 1
	W = [-ω^2             0          -im*ω*Pe*ϕ0;
		   -ω^2*ϕa*DD      -ω^2*DS     -im*ω*Pe*(ϕa*s+DS);
		   -im*ω*Pe*ϕa*DSP -im*ω*Pe*DS -ω^2*DS-2]
	values, vectors = eigen(W)
	return real(values[3])
end

function relambw2(w; ϕa = 0.1, ϕ = 0.7, v = 7.5)
	ω = w;
	Pe = v;
	ϕ0 = 1 - ϕ;
	DS = self_diff(ϕ);
	DSP = self_diff_prime(ϕ);

	DD = (1-DS)/ϕ
	s = DD - 1
	W = [-ω^2             0          -im*ω*Pe*ϕ0;
		   -ω^2*ϕa*DD      -ω^2*DS     -im*ω*Pe*(ϕa*s+DS);
		   -im*ω*Pe*ϕa*DSP -im*ω*Pe*DS -ω^2*DS-2]
	values, vectors = eigen(W)
	return real(values[2])
end

function relambw3(w; ϕa = 0.1, ϕ = 0.7, v = 7.5)
	ω = w;
	Pe = v;
	ϕ0 = 1 - ϕ;
	DS = self_diff(ϕ);
	DSP = self_diff_prime(ϕ);

	DD = (1-DS)/ϕ
	s = DD - 1
	W = [-ω^2             0          -im*ω*Pe*ϕ0;
		   -ω^2*ϕa*DD      -ω^2*DS     -im*ω*Pe*(ϕa*s+DS);
		   -im*ω*Pe*ϕa*DSP -im*ω*Pe*DS -ω^2*DS-2]
	values, vectors = eigen(W)
	return real(values[1])
end

function imlambw(w; ϕa = 0.1, ϕ = 0.7, v = 7.5)
	ω = w;
	Pe = v;
	ϕ0 = 1 - ϕ;
	DS = self_diff(ϕ);
	DSP = self_diff_prime(ϕ);

	DD = (1-DS)/ϕ
	s = DD - 1
	W = [-ω^2             0          -im*ω*Pe*ϕ0;
		   -ω^2*ϕa*DD      -ω^2*DS     -im*ω*Pe*(ϕa*s+DS);
		   -im*ω*Pe*ϕa*DSP -im*ω*Pe*DS -ω^2*DS-2]
	values, vectors = eigen(W)
	values = sort(values, by = x->imag.(x))
	return imag(values[3])
end

function imlambw2(w; ϕa = 0.1, ϕ = 0.7, v = 7.5)
	ω = w;
	Pe = v;
	ϕ0 = 1 - ϕ;
	DS = self_diff(ϕ);
	DSP = self_diff_prime(ϕ);

	DD = (1-DS)/ϕ
	s = DD - 1
	W = [-ω^2             0          -im*ω*Pe*ϕ0;
		   -ω^2*ϕa*DD      -ω^2*DS     -im*ω*Pe*(ϕa*s+DS);
		   -im*ω*Pe*ϕa*DSP -im*ω*Pe*DS -ω^2*DS-2]
	values, vectors = eigen(W)
	values = sort(values, by = x->imag.(x))
	return imag(values[1])
end

function imlambw3(w; ϕa = 0.1, ϕ = 0.7, v = 7.5)
	ω = w;
	Pe = v;
	ϕ0 = 1 - ϕ;
	DS = self_diff(ϕ);
	DSP = self_diff_prime(ϕ);

	DD = (1-DS)/ϕ
	s = DD - 1
	W = [-ω^2             0          -im*ω*Pe*ϕ0;
		   -ω^2*ϕa*DD      -ω^2*DS     -im*ω*Pe*(ϕa*s+DS);
		   -im*ω*Pe*ϕa*DSP -im*ω*Pe*DS -ω^2*DS-2]
	values, vectors = eigen(W)
	values = sort(values, by = x->imag.(x))
	return imag(values[2])
end


function is_stable_value(ϕa, ϕp; Pe = 10.0)

	ϕ = ϕa + ϕp
	max_val = relambw(0.01; ϕa = ϕa, ϕ = ϕ, v = Pe)

	# get extreme values 
	f(x) = Dw(ϕ, ϕa, Pe, x, ds(ϕ), dsp(ϕ))*Dw2(ϕ, ϕa, Pe, x, ds(ϕ), dsp(ϕ))*Dw3(ϕ, ϕa, Pe, x, ds(ϕ), dsp(ϕ))
	try
		ws = find_zeros(f, (1e-10, 10))
		for w in ws
			max_val = max(max_val, relambw(w; ϕa = ϕa, ϕ = ϕ, v = Pe))
		end
	catch
	end

	return max_val-1e-10
end

function is_stable_value_real(ϕa, ϕp; Pe = 10.0)
	ϕ = ϕa + ϕp
	expr2 = (β(ϕa, ϕp) + α(ϕa, ϕp)*Pe^2)
	expr1 = β(ϕa, ϕp)-α(ϕa, ϕp)*Pe^2-4*ϕ
	return real(expr1+sqrt(complex(expr2^2 - 4*ϕp*β(ϕa, ϕp)*ds(ϕ)^2*Pe^2)))
end

function is_stable_value_finite(ϕa, ϕp; Pe = 10.0, Lx = 20)
	ω = 2*π/Lx;
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
	return real(values[3])
end

function is_complex_value(ϕa, ϕp; Pe = 10.0) # positive implies complex
	ϕ = ϕa + ϕp
	# get extreme values 
	f(x) = Dw(ϕ, ϕa, Pe, x, ds(ϕ), dsp(ϕ))*Dw2(ϕ, ϕa, Pe, x, ds(ϕ), dsp(ϕ))*Dw3(ϕ, ϕa, Pe, x, ds(ϕ), dsp(ϕ))
	ws = find_zeros(f, (1e-10, 10))
	try
		if ws == [1e-10]
			expr2 = (β(ϕa, ϕp) + α(ϕa, ϕp)*Pe^2)
			return -(expr2^2 - 4*ϕp*β(ϕa, ϕp)*ds(ϕ)^2*Pe^2)
		else
			ind = argmax(relambw.(ws; ϕa = ϕa, ϕ = ϕ, v = Pe))
			return discriminant_f(ϕ, ϕa, Pe, ws[ind], ds(ϕ), dsp(ϕ)) - 1e-15
		end
	catch
		# expr2 = (β(ϕa, ϕp) + α(ϕa, ϕp)*Pe^2)
		return -1e-22 #-(expr2^2 - 4*ϕp*β(ϕa, ϕp)*ds(ϕ)^2*Pe^2)
	end
end

function return_binodal(; Pe = 7.5)
	filename = datadir("binodal", "Pe=$(Pe).jld2")

	try
		data = wload(filename)
		println("loaded: ", "Pe=$(Pe)")
		@unpack Pe, γs, ϕ1s, ϕ2s = data
		return γs, ϕ1s, ϕ2s
	catch
		println("computing binodal: ", "Pe=$(Pe)")

		ϕ = find_zero(x->(pspin.(x; Pe = Pe)/(1-x) - pspin.(x+1e-10; Pe = Pe)/(1-x-1e-10))/1e-10, (0.6, 0.99))
		γ_max = (1-ϕ+pspin.(ϕ; Pe = Pe))/(1-ϕ)

		γs, ϕ1s, ϕ2s = [], [], []

		for γ in 1:0.001:(γ_max-0.001)
			try
				find_sol, lower_limits, upper_limits = colapse_sol_interval(; Pe = Pe, γ = γ, initial_Δ = 1e-6, max_iter = 40, tol = 1e-8, atol = 1e-12)
				if find_sol
					ϕ1, ϕ2 = lower_limits[1], upper_limits[2]
					push!(ϕ1s, ϕ1)
					push!(ϕ2s, ϕ2)
					push!(γs, γ)
				end
			catch
				break
			end
		end

		for γ in (γ_max-0.001):0.0001:(γ_max-0.0001)
			try
				find_sol, lower_limits, upper_limits = colapse_sol_interval(; Pe = Pe, γ = γ, initial_Δ = 1e-8, max_iter = 200, tol = 1e-8, atol = 1e-12)
				if find_sol
					ϕ1, ϕ2 = lower_limits[1], upper_limits[2]
					push!(ϕ1s, ϕ1)
					push!(ϕ2s, ϕ2)
					push!(γs, γ)
				end
			catch
				break
			end
		end

		for γ in (γ_max-0.0001):0.00001:(γ_max-0.00001)
			try
				find_sol, lower_limits, upper_limits = colapse_sol_interval(; Pe = Pe, γ = γ, initial_Δ = 1e-8, max_iter = 200, tol = 1e-8, atol = 1e-12)
				if find_sol
					ϕ1, ϕ2 = lower_limits[1], upper_limits[2]
					push!(ϕ1s, ϕ1)
					push!(ϕ2s, ϕ2)
					push!(γs, γ)
				end
			catch
				break
			end
		end

		for γ in (γ_max-0.00001):0.000001:(γ_max)
			try
				find_sol, lower_limits, upper_limits = colapse_sol_interval(; Pe = Pe, γ = γ, initial_Δ = 1e-8, max_iter = 200, tol = 1e-8, atol = 1e-12)
				if find_sol
					ϕ1, ϕ2 = lower_limits[1], upper_limits[2]
					push!(ϕ1s, ϕ1)
					push!(ϕ2s, ϕ2)
					push!(γs, γ)
				end
			catch
				break
			end
		end

		data = Dict{String, Any}()
		@pack! data = Pe, γs, ϕ1s, ϕ2s
		safesave(filename, data)
		println("saved: ", "Pe=$(Pe)")
		return γs, ϕ1s, ϕ2s

	end
end
function return_spin(; Pe = Pe, Δϕ = Δϕ, ϕp_max = 1, save_on = true)
	try
		file_name = datadir("spinodal", "Δϕ=$(Δϕ)_Pe=$(Pe).jld2")
		data = load(file_name)
		@unpack ϕas_left, ϕas_right, ϕps, indl, indr = data
		nϕps = [z for (x, y, z) in zip(ϕas_left, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		nϕas_right = [y for (x, y, z) in zip(ϕas_left, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		nϕas_left = [x for (x, y, z) in zip(ϕas_left, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		return nϕas_left, nϕas_right, nϕps, indl, indr
	catch
		println("computing: ", "Δϕ=$(Δϕ)_Pe=$(Pe)")
		# Sample the spinodal horizontally
		# ϕp_grid = Δϕ:Δϕ:ϕp_max
		ϕp_grid = [0.001; Δϕ:Δϕ:ϕp_max]
		ϕps = []
		ϕas_left = []
		ϕas_right = []
		for ϕp in ϕp_grid
			println("doing ϕp = ", ϕp)
			f(x) = is_stable_value(x, ϕp; Pe = Pe)
			try
				# ϕal, ϕar = find_zeros(f,(1e-10,1-ϕp-1e-10))
				ϕal, ϕar = find_zeros(f, (0.3, 0.9))
				println("found [ϕal, ϕar, ϕp] = [", d3(ϕal), ", ", d3(ϕar), ", ", d3(ϕp), "]")
				push!(ϕas_left, ϕal)
				push!(ϕas_right, ϕar)
				push!(ϕps, ϕp)
			catch
				println("no zeros for ϕp = ", ϕp)
			end
		end
		# Sample the spinodal vertically to resolve the tip of the spinodal
		if ϕas_left == [] && ϕas_right == []
			println("no spinodal found for this grid")
			return [], [], [], 0, 0  # Return empty arrays and default indices
		else
			ϕa_grid = ϕas_left[end]+0.001:Δϕ:ϕas_right[end]-0.001
			for ϕa in ϕa_grid
				println("doing ϕa = ", d3(ϕa))
				g(y) = is_stable_value(ϕa, y; Pe = Pe)
				try
					ϕp = find_zero(g, (1e-10, 1-ϕa-1e-8))
					println("found [ϕa, ϕp] = [", d3(ϕa), ", ", d3(ϕp), "]")
					f(x) = is_stable_value(x, ϕp; Pe = Pe)
					ϕal, ϕar = find_zeros(f, (0, 1-ϕp-1e-8))
					println("found [ϕal, ϕar, ϕp] = [", d3(ϕal), ", ", d3(ϕar), ", ", d3(ϕp), "]")
					push!(ϕas_left, ϕal)
					push!(ϕas_right, ϕar)
					push!(ϕps, ϕp)
				catch
					println("no zeros for ϕa = ", d3(ϕa))
				end
			end

			# Find the points where the nature of the instability changes (from real to complex)
			indl = argmax(sign.(is_complex_value.(ϕas_left, ϕps; Pe = Pe)))
			indr = argmax(sign.(is_complex_value.(ϕas_right, ϕps; Pe = Pe)))

			file_name = datadir("spinodal", "Δϕ=$(Δϕ)_Pe=$(Pe).jld2")
			data = Dict{String, Any}()
			@pack! data = ϕas_left, ϕas_right, ϕps, indl, indr
			safesave(file_name, data)
			println("saved: ", "Δϕ=$(Δϕ)_Pe=$(Pe)")
			return ϕas_left, ϕas_right, ϕps, indl, indr
		end
	end
end

function return_spin_finite(; Pe = 7.5, Δϕ = 0.01, Lx = 20)
	ϕp_grid = Δϕ:Δϕ:1
	ϕps = []
	ϕas_left = []
	ϕas_right = []
	for ϕp in ϕp_grid
		f(x) = is_stable_value_finite(x, ϕp; Pe = Pe, Lx = Lx)
		try
			ϕal, ϕar = find_zeros(f, (0, 1-ϕp-1e-8))
			push!(ϕas_left, ϕal)
			push!(ϕas_right, ϕar)
			push!(ϕps, ϕp)
		catch
		end
	end
	return ϕas_left, ϕas_right, ϕps
end

function return_spin_from_grid(grid_name; max_ϕa = 1.0, Pe = Pe, γ_grid = 0.1:0.1:1, ϕ1_grid = [], ϕ2_grid = [], ϕp_grid = 0.1:0.1:1)
	try
		file_name = datadir("spinodal", "$(grid_name).jld2")
		data = load(file_name)
		@unpack ϕas_left, ϕas_right, ϕps, γs, ϕ1s, ϕ2s = data
		nϕps = [z for (x, y, z) in zip(ϕas_left, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		nϕas_right = [y for (x, y, z) in zip(ϕas_left, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		nϕas_left = [x for (x, y, z) in zip(ϕas_left, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		nϕas_left = [x for (x, y, z) in zip(ϕas_left, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		γs = [x for (x, y, z) in zip(γs, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		ϕ1s = [x for (x, y, z) in zip(ϕ1s, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]
		ϕ2s = [x for (x, y, z) in zip(ϕ2s, ϕas_right, ϕps) if (y>0.5)|(z > 0.2)]

		return nϕas_left, nϕas_right, nϕps, γs, ϕ1s, ϕ2s
	catch
		println("computing: ", "$(grid_name)")
		ϕps = []
		ϕas_left = []
		ϕas_right = []
		γs = []
		ϕ1s = []
		ϕ2s = []
		for (ϕp, γ, ϕ1, ϕ2) in zip(ϕp_grid, γ_grid, ϕ1_grid, ϕ2_grid)
			f(x) = is_stable_value(x, ϕp; Pe = Pe)
			try
				# ϕal, ϕar = find_zeros(f,(0,1-ϕp-1e-8))
				ϕal, ϕar = find_zeros(f, (0.3, 0.9))
				push!(ϕas_left, ϕal)
				push!(ϕas_right, ϕar)
				push!(ϕps, ϕp)
				push!(γs, γ)
				push!(ϕ1s, ϕ1)
				push!(ϕ2s, ϕ2)
			catch

				push!(ϕas_left, max_ϕa)
				push!(ϕas_right, max_ϕa)
				push!(ϕps, ϕp)
				push!(γs, γ)
				push!(ϕ1s, ϕ1)
				push!(ϕ2s, ϕ2)
			end
		end

		file_name = datadir("spinodal", "$(grid_name).jld2")
		data = Dict{String, Any}()
		@pack! data = ϕas_left, ϕas_right, ϕps, γs, ϕ1s, ϕ2s
		safesave(file_name, data)
		println("saved: ", "$(grid_name)")
		return ϕas_left, ϕas_right, ϕps, γs, ϕ1s, ϕ2s
	end
end


function return_spin_from_grid_real(; max_ϕa = 1.0, Pe = Pe, γ_grid = 0.1:0.1:1, ϕ1_grid = [], ϕ2_grid = [], ϕp_grid = 0.1:0.1:1)
	ϕps = []
	ϕas_left = []
	ϕas_right = []
	γs = []
	ϕ1s = []
	ϕ2s = []
	for (ϕp, γ, ϕ1, ϕ2) in zip(ϕp_grid, γ_grid, ϕ1_grid, ϕ2_grid)
		f(x) = is_stable_value_real(x, ϕp; Pe = Pe)
		try
			ϕal, ϕar = find_zeros(f, (0, 1-ϕp-1e-8))
			# ϕal, ϕar = find_zeros(f,(0.3,0.9))
			push!(ϕas_left, ϕal)
			push!(ϕas_right, ϕar)
			push!(ϕps, ϕp)
			push!(γs, γ)
			push!(ϕ1s, ϕ1)
			push!(ϕ2s, ϕ2)
		catch

			push!(ϕas_left, max_ϕa)
			push!(ϕas_right, max_ϕa)
			push!(ϕps, ϕp)
			push!(γs, γ)
			push!(ϕ1s, ϕ1)
			push!(ϕ2s, ϕ2)
		end
	end

	return ϕas_left, ϕas_right, ϕps, γs, ϕ1s, ϕ2s
end
#

function return_complex_boundary_pt(ϕa; Pe = 10.0)
	f(x) = is_complex_value(ϕa, x; Pe = Pe)
	return find_zeros(f, (0, 1-ϕa-1e-8))
end

function return_complex_boundary_outer(ϕas; Pe = 10.0)
	high_ϕps = []
	low_ϕps = []
	ϕa_sols = []
	for ϕa in ϕas
		ϕp = return_complex_boundary_pt(ϕa; Pe = Pe)
		if length(ϕp)==2
			push!(high_ϕps, ϕp[2])
			push!(low_ϕps, ϕp[1])
			push!(ϕa_sols, ϕa)
		end
	end
	return ϕa_sols, low_ϕps, high_ϕps
end

function return_complex_boundary_inner(ϕas; Pe = 10.0)
	high_ϕps1 = []
	high_ϕps2 = []
	low_ϕps1 = []
	low_ϕps2 = []
	ϕa_sols = []
	for ϕa in ϕas
		ϕp = return_complex_boundary_pt(ϕa; Pe = Pe)
		if length(ϕp)==4
			push!(high_ϕps1, ϕp[2])
			push!(low_ϕps1, ϕp[1])
			push!(high_ϕps2, ϕp[4])
			push!(low_ϕps2, ϕp[3])
			push!(ϕa_sols, ϕa)
		end
	end
	return ϕa_sols, low_ϕps1, high_ϕps1, low_ϕps2, high_ϕps2
end

#

function return_stable_boundary_pt(ϕa; Pe = 10.0)
	f(x) = is_stable_value(ϕa, x; Pe = Pe)
	return find_zeros(f, (0, 1-ϕa-1e-8))
end

function return_stable_boundary_outer(ϕas; Pe = 10.0)
	high_ϕps = []
	low_ϕps = []
	ϕa_low = []
	ϕa_high = []
	for ϕa in ϕas
		ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
		if length(ϕp)>1
			if (is_complex_value(ϕa, minimum(ϕp); Pe = Pe)<0)
				push!(low_ϕps, minimum(ϕp))
				push!(ϕa_low, ϕa)
			end
			if (is_complex_value(ϕa, maximum(ϕp); Pe = Pe)<0)
				push!(high_ϕps, maximum(ϕp))
				push!(ϕa_high, ϕa)
			end
		elseif length(ϕp)>0
			if (is_complex_value(ϕa, ϕp[1]; Pe = Pe)<0)
				push!(high_ϕps, ϕp[1])
				push!(ϕa_high, ϕa)
			end
		end
	end
	return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary(ϕas; Pe = 10.0)
	high_ϕps = []
	low_ϕps = []
	ϕa_low = []
	ϕa_high = []
	for ϕa in ϕas
		ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
		if length(ϕp)==2
			push!(low_ϕps, minimum(ϕp))
			push!(ϕa_low, ϕa)
			push!(high_ϕps, maximum(ϕp))
			push!(ϕa_high, ϕa)
		elseif length(ϕp)==1
			push!(high_ϕps, ϕp[1])
			push!(ϕa_high, ϕa)
		end
	end
	return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary_outer(ϕas; Pe = 10.0)
	high_ϕps = []
	low_ϕps = []
	ϕa_low = []
	ϕa_high = []
	for ϕa in ϕas
		ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
		if length(ϕp)>1
			if (is_complex_value(ϕa, minimum(ϕp); Pe = Pe)<0)
				push!(low_ϕps, minimum(ϕp))
				push!(ϕa_low, ϕa)
			end
			if (is_complex_value(ϕa, maximum(ϕp); Pe = Pe)<0)
				push!(high_ϕps, maximum(ϕp))
				push!(ϕa_high, ϕa)
			end
		elseif length(ϕp)>0
			if (is_complex_value(ϕa, ϕp[1]; Pe = Pe)<0)
				push!(high_ϕps, ϕp[1])
				push!(ϕa_high, ϕa)
			end
		end
	end
	return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary_inner(ϕas; Pe = 10.0)
	high_ϕps = []
	low_ϕps = []
	ϕa_low = []
	ϕa_high = []
	for ϕa in ϕas
		ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
		if length(ϕp)==2
			#if (is_complex_value(ϕa, minimum(ϕp); Pe = Pe)>0)
			push!(low_ϕps, minimum(ϕp))
			push!(ϕa_low, ϕa)
			#end
			#if (is_complex_value(ϕa, maximum(ϕp); Pe = Pe)>0)
			push!(high_ϕps, maximum(ϕp))
			push!(ϕa_high, ϕa)
			#end
		elseif length(ϕp)>0
			#if (is_complex_value(ϕa, ϕp[1]; Pe = Pe)>0)
			push!(high_ϕps, ϕp[1])
			push!(ϕa_high, ϕa)
			#end
		end
	end
	return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary_extra(ϕas; Pe = 10.0)
	high_ϕps = []
	mid_ϕps = []
	low_ϕps = []
	ϕa_low = []
	ϕa_mid = []
	ϕa_high = []
	for ϕa in ϕas
		ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
		if length(ϕp)==3
			push!(low_ϕps, minimum(ϕp))
			push!(ϕa_low, ϕa)
			push!(high_ϕps, maximum(ϕp))
			push!(ϕa_high, ϕa)
			push!(mid_ϕps, ϕp[2])
			push!(ϕa_mid, ϕa)
		end
	end
	return ϕa_low, ϕa_mid, ϕa_high, low_ϕps, mid_ϕps, high_ϕps
end


function is_fake_value(ϕa, ϕp; Pe = 10.0)
	ϕ  = ϕa + ϕp
	ϕ0 = 1 - ϕ
	ds  = self_diff(ϕ)
	dsp = self_diff_prime(ϕ)
	return Pe^2 + 2/(ds*ϕp + ϕ0*(ds + dsp*ϕa))
end

function fake_spin_pt(ϕa; Pe = 10.0)
	f(x) = is_fake_value(ϕa, x; Pe = Pe)
	return find_zeros(f, (0.0, 1-ϕa-1e-8))[1]
end

function fake_spin_boundary(ϕas; Pe = 10.0)
	ϕp_pts = []
	ϕa_pts = []
	for ϕa in ϕas
		try
			ϕp = fake_spin_pt(ϕa; Pe = Pe)
			push!(ϕa_pts, ϕa)
			push!(ϕp_pts, ϕp)
		catch
		end
	end
	return ϕa_pts, ϕp_pts
end

#plot approx spinodal
function pm_lin_real_value(ϕa, ϕp; Pe = 7.5, Dθ = 400.0)
	ω = 2*π/sqrt(Dθ);
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
	return real(values[3])
end

function return_finite_stable_boundary_pt(ϕa; Pe = 10.0, Dθ = 400.0)
	f(x) = pm_lin_real_value(ϕa, x; Pe = Pe, Dθ = Dθ)
	return find_zeros(f, (0, 1-ϕa-1e-8))
end

function return_finite_stable_boundary(ϕas; Pe = 10.0, Dθ = 400.0)
	high_ϕps = []
	low_ϕps = []
	ϕa_low = []
	ϕa_high = []
	for ϕa in ϕas
		ϕp = return_finite_stable_boundary_pt(ϕa; Pe = Pe, Dθ = Dθ)
		if length(ϕp)==2
			push!(low_ϕps, minimum(ϕp))
			push!(ϕa_low, ϕa)
			push!(high_ϕps, maximum(ϕp))
			push!(ϕa_high, ϕa)
		elseif length(ϕp)==1
			push!(high_ϕps, ϕp[1])
			push!(ϕa_high, ϕa)
		end
	end
	return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_finite_stable_boundary_extra(ϕas; Pe = 10.0, Dθ = 400.0)
	high_ϕps = []
	mid_ϕps = []
	low_ϕps = []
	ϕa_low = []
	ϕa_mid = []
	ϕa_high = []
	for ϕa in ϕas
		ϕp = return_finite_stable_boundary_pt(ϕa; Pe = Pe, Dθ = Dθ)
		if length(ϕp)==3
			push!(low_ϕps, minimum(ϕp))
			push!(ϕa_low, ϕa)
			push!(high_ϕps, maximum(ϕp))
			push!(ϕa_high, ϕa)
			push!(mid_ϕps, ϕp[2])
			push!(ϕa_mid, ϕa)
		end
	end
	return ϕa_low, ϕa_mid, ϕa_high, low_ϕps, mid_ϕps, high_ϕps
end

#
