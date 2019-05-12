#=
main:
- Julia version:
- Author: quangio
- Date: 2019-05-06
=#



using JuMP
using SCIP
using DelimitedFiles
include("input.jl")

model = Model(with_optimizer(SCIP.Optimizer))


dist_mat_gr17 = [
  0 633 257  91 412 150  80 134 259 505 353 324  70 211 268 246 121
633   0 390 661 227 488 572 530 555 289 282 638 567 466 420 745 518
257 390   0 228 169 112 196 154 372 262 110 437 191  74  53 472 142
 91 661 228   0 383 120  77 105 175 476 324 240  27 182 239 237  84
412 227 169 383   0 267 351 309 338 196  61 421 346 243 199 528 297
150 488 112 120 267   0  63  34 264 360 208 329  83 105 123 364  35
 80 572 196  77 351  63   0  29 232 444 292 297  47 150 207 332  29
134 530 154 105 309  34  29   0 249 402 250 314  68 108 165 349  36
259 555 372 175 338 264 232 249   0 495 352  95 189 326 383 202 236
505 289 262 476 196 360 444 402 495   0 154 578 439 336 240 685 390
353 282 110 324  61 208 292 250 352 154   0 435 287 184 140 542 238
324 638 437 240 421 329 297 314  95 578 435   0 254 391 448 157 301
 70 567 191  27 346  83  47  68 189 439 287 254   0 145 202 289  55
211 466  74 182 243 105 150 108 326 336 184 391 145   0  57 426  96
268 420  53 239 199 123 207 165 383 240 140 448 202  57   0 483 153
246 745 472 237 528 364 332 349 202 685 542 157 289 426 483   0 336
121 518 142  84 297  35  29  36 236 390 238 301  55  96 153 336   0] # shoud return 2085
dist_mat_5 = [
0.0  3.0  4.0  2.0  7.0
3.0  0.0  4.0  6.0  3.0
4.0  4.0  0.0  5.0  8.0
2.0  6.0  5.0  0.0  6.0
7.0  3.0  8.0  6.0  0.0
] # shoud return 19.0

n = 22 # number of vertices
c = TestSet.eil22_capacities / 100  # cost assigned to each vertex (weight needed to be picked up)
d = TestSet.loadVrpEil22() # distance matrix of the graph

w0 = 100
no_salesman = 5
cap = 300

range = collect(1:n)

for i = 1:n
    d[i, i] = 0
    # c[i] = rand() * 100
end
c[1] = 0

println(string("solving TSP for capacities: ", c))

@variable(model, y[1:n, 1:n])
@variable(model, x[1:n, 1:n], Bin)

@objective(model, Min, sum(y[i, j] * d[i, j] for i = 1:n, j = range[1:n .!= i]))

@constraint(model, sum(x[1, i] for i = 2:n) == no_salesman)
@constraint(model, sum(x[i, 1] for i = 2:n) == no_salesman)

for i = 2:n
    @constraint(model, y[1, i] == w0 * x[1, i])
    @constraint(model, sum(x[i, j] for j = range[1:n .!= i]) == 1)
    @constraint(model, sum(x[j, i] for j = range[1:n .!= i]) == 1)
    @constraint(model, c[i] == sum(y[i, j] - y[j, i] for j = range[1:n .!= i]))
end

for i = 1:n
    for j = 1:n
        if i == j continue end
        @constraint(model, y[i, j] >= (w0 + c[i]) * x[i, j])
        @constraint(model, y[i, j] <= (w0 + cap - c[j]) * x[i, j])
    end
end


function solved(m)
    x_val = JuMP.value.(x)
    cycle_idx = [1] # we can adapt this to support mTSP, just keep this for simplicity for now
    while true
        v, idx = findmax(x_val[cycle_idx[length(cycle_idx)], :])
        if idx == 1
            break
        else
            push!(cycle_idx, idx)
        end
    end
    if length(cycle_idx) < n
        @constraint(m, sum(x[i, j] for i = cycle_idx, j = cycle_idx[cycle_idx .!= i]) <= length(cycle_idx) - 1)
        return false
    end
    return true
end

if no_salesman == 0
    JuMP.optimize!(model)
    while !solved(model)
        JuMP.optimize!(model)
    end
end

if no_salesman >= 1000
    @variable(model, 0 <= u[1:n] <= n - 1)
    for i = 2:n
        for j = collect(2:n)[2:n .!= i]
            @constraint(model, u[i] - u[j] + n * x[i, j] <= n - 1)
        end
    end
    @time optimize!(model)
end

optimize!(model)

println(string("Objective value:", JuMP.objective_value(model)))

println("Arcs: ")

for i = range
    for j = range
        if value(x[i, j]) > 0.0
            println(string(i, "->", j, "[label=", value(y[i, j]), "]")) # you can use graphviz to visualize this
        end
    end
end

println("Arcs with capacity: ")

for i = range
    for j = range
        if value(x[i, j]) > 0.0
            println(string(floor(c[i]), "->", floor(c[j]), "[label=", floor(value(y[i, j])), "]")) # you can use graphviz to visualize this
        end
    end
end

# open("model2.lp", "w") do f
#    print(f, model)
# end
