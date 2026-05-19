using Roots

include("startup.jl")


M = 20
grid = [[i/M for i in 0:M] for _ in 1:1]
cov = 0.2#0.2 #0.1
alpha = 0.8 


# Transition dynmics
function reflect(x)
	count = 0
	while (x < 0) || (x > 1)
		if x < 0
			x = -1 * x
		end
		if x > 1
			x = 2 - x
		end
		count += 1 
	end
	if count > 1
		println("reflection count = $count ")
	end
	return x
end

function sample_next(x)
	r = Normal(0, 1)
	y = x[1] + alpha*(0.5-x[1]) + cov*rand(r, 1)[1]
	return reflect(y)
end

lowerbound = pdf(alpha*0.5 + cov*Normal(0,1), 1.0) + pdf(alpha*0.5 + cov*Normal(0,1), -1.0)
println("The lower bound gamma = $lowerbound")


#=
# These boundary conditions dont result in a Gaussian
#
function sample_next(x)
	r = Normal(0, 1)
	y = x[1] + alpha*(0.5-x[1]) + cov*rand(r, 1)[1]
	y = min(y, 1.0)
	y = max(y, 0.0)
	return y
end

=#
# Initial Condition
x = 0.7

# Run it for data
println("creating data")
data = [[x]]
D = 6 
for t in 1:(10^D)
	push!(data, [sample_next(data[end])])
end

# If you want the histogram version of the data;
do_histogram = true
if do_histogram
	println("plotting histogram")
	datah = [d[1] for d in data]
	ph = Plots.histogram(datah, normalize= :pdf)
	display(ph)
end

# Here is the Theoretical Stationary Distribtution
gram = cov^2 * (1 / (1-(1-alpha)^2)) 
plotpdf = Normal(0.5, sqrt(gram))
renorm_fac = cdf(plotpdf,1.0) - cdf(plotpdf, 0.0)

function thr_pdf(x) 
	return pdf(plotpdf, x) + pdf(plotpdf, -1*x) + pdf(plotpdf, 2-x) 
end

function pdf_derivative(x)
	return (1 / sqrt(2*pi*gram))*(-2*(x-0.5)/(2*gram))*exp((-1*(x-0.5)^2)/(2*gram))
end

function pdf_derivative2(x)
	return (1 / sqrt(2*pi*gram))*((((x-0.5)^2)/(gram^2))-(1/gram))*exp((-1*(x-0.5)^2)/(2*gram))
end

function reflected_derv(x)
	return pdf_derivative(x) - pdf_derivative(-1*x) - pdf_derivative(2-x)
end

function reflected_derv2(x)
	return pdf_derivative2(x) + pdf_derivative2(-1*x) +  pdf_derivative2(2-x)
end





Plots.plot!(thr_pdf, 0.0:0.001:1.0)
display(ph)


# Derivative nonsense
#=
pd = Plots.plot(reflected_derv, 0:0.01:1)
rt1 = find_zero(reflected_derv2, 0.5 - sqrt(gram))
rt2 = find_zero(reflected_derv2, 0.5 + sqrt(gram))
println("rts = $([rt1, rt2])")
Plots.plot!([rt1], [reflected_derv(rt1)], markershape=:circle)
Plots.plot!([rt2], [reflected_derv(rt2)], markershape=:circle)
println("rts = $([rt1, rt2])")
display(pd)
=#

# Here is the function for the kernel density
function kern_dens(x, y)
	# this is a gaussian reflected at the borders 0, 1
	random_var = x + alpha*(0.5 - x) + cov*Normal(0,1) 
	return pdf(random_var, y) + pdf(random_var, -1*y) + pdf(random_var, 2-y)
end

#finegrid = [i for i in 0:0.01:1]
#kern = [kern_dens(x, y) for x in finegrid, y in finegrid]
#pheat_true = Plots.heatmap(finegrid, finegrid, kern)



println("type 1 to continue. Else exit")
#r = readline()
if true#r == "1"

	# Frankliin Approximation!

	# 1D basis for the stationary distribution
	println("1D basis approximation")
	numbasis = 2^4 + 1#2^6
	basis = @time Func.create_basis(1, numbasis)
	println("time for emp app")
	em1_c = @time Func.approximate_emp(basis, data) 
	em1 = em1_c[1]
	cem1 = em1_c[2]

	nem1 = Func.normalize_1(em1) 

	# plot the approximation of the stationary
	pn = Plots.plot!(x -> Func.eval_f(nem1, [x]), 0:0.01:1, label="Approximation")
	display(pn)
end
println("type 1 to continue. Else exit")
#r = readline()
if true#r == "1"
	# Use the stationary approximation for the kernel approximation

	# construct the paired dataset
	data2 = [[data[i][1], data[i+1][1]] for i in 1:length(data)-1]
	println("time for 2d basis")
	basis2 = @time Func.create_basis(2, 8)#16)

	println("time for kern app")
	
	em2_c = @time Func.approximate_kern(basis2, data2, em1)
	em2 = em2_c[1]
	cem2 = em2_c[2]

	em2 = Func.normalize_1(em2)
	
	z = Func.plot(em2)
	Makie.surface(z...)
	println(minimum(z[3]))

end






