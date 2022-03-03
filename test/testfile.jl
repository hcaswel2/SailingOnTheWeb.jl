using SailingOnTheWeb

β1_1_lst = [55,125,-125,-55] .* pi/180
α1_1_lst = [0,30,30,0] .* pi/180
u_o1 = [0.0; -55]  

β1_2_lst = [62,118,-122,-62] .* pi/180
α1_2_lst = [10.32,52,52,0.3] .* pi/180
u_o2 = [0.0; -61.64] 

# Define Constants
h, m, h1, h2, Cᴰ, Cᴸ, r_island, ut_o, t_int, A, h_avg, ρ = def_constants()

# Calculate/generate wind
VWx_avg, VWy_avg, VW_avg = generate_wind(h1,h2,h_avg,r_island)

# Package Constants
package_1 = [α1_1_lst[1], α1_1_lst, β1_1_lst[1], β1_1_lst, 1, VWx_avg, VWy_avg, m, h, A, Cᴰ, Cᴸ, 0, 0,r_island];
package_2 = [α1_2_lst[1], α1_2_lst, β1_2_lst[1], β1_2_lst, 1, VWx_avg, VWy_avg, m, h, A, Cᴰ, Cᴸ, 0, 0,r_island];

# Define Eq. of motion


# Define callbacks
cbs1, cbs2 = def_callbacks()

# Solve Equation of motion for boat 1
problem1 = SecondOrderODEProblem(eq_motion, ut_o, u_o1, (t_int[1],t_int[end]),package_1)
solution1 = solve(problem1, DPRKN8() ,tstops=t_int,callback=cbs1,save_everystep=true)
u1 = [(x->x[1]).(solution1.u), (x->x[2]).(solution1.u),(x->x[3]).(solution1.u),(x->x[4]).(solution1.u)];
t_1 = solution1.t
t_β1 = saved_values1.t #retrieve tracked boat angle timestep
α1_1list = (x->x[1]).(saved_values1.saveval) #retrieve tracked sail angle
β1_1list = (x->x[2]).(saved_values1.saveval) #retrieve tracked boat angle

# Solve Equation of motion for boat 2
problem2 = SecondOrderODEProblem(eq_motion, ut_o, u_o2, (t_int[1],t_int[end]),package_2,save_everystep=true)
solution2 = solve(problem2, DPRKN8(),tstops=t_int,callback=cbs2)
u2 = [(x->x[1]).(solution2.u), (x->x[2]).(solution2.u),(x->x[3]).(solution2.u),(x->x[4]).(solution2.u)];
t_2 = solution2.t 
t_β2 = saved_values2.t #retrieve tracked boat angle timestep
α1_2list = (x->x[1]).(saved_values2.saveval) #retrieve tracked sail angle
β1_2list = (x->x[2]).(saved_values2.saveval) #retrieve tracked boat angle

# Animate

anim_boat(u1,β1_1list,α1_1list,u2,β1_2list,α1_2list,VWx_avg,VWy_avg,VW_avg,U_avg,t_1,t_2,t_β1,t_β2,5)