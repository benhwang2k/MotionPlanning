numbasis = 2^4 + 1#2^6
basis = @time Func.create_basis(1, numbasis)
for b in basis
  # plot the basis 
  Plots.plot!(x -> Func.eval_f(b, [x]), 0:0.01:1, label="Approximation")
end

display(pn)
