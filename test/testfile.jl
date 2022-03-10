using SailingOnTheWeb

β1_1_lst = [55,125,-125,-55] .* pi/180
α1_1_lst = [0,30,30,0] .* pi/180
u_o1 = [0.0; -55]  

β1_2_lst = [62,118,-122,-62] .* pi/180
α1_2_lst = [10.32,52,52,0.3] .* pi/180
u_o2 = [0.0; -61.64] 

t_int = 0.0:0.1:1

# Define Constants
w, h, m, h1, h2, Cᴰ, Cᴸ, r_island, ut_o, A, h_avg, ρ = def_constants()

# Calculate/generate wind
VWx_avg, VWy_avg, VW_avg = generate_wind(h1,h2,h_avg,r_island)

# Package Constants
package_1 = [α1_1_lst[1], α1_1_lst, β1_1_lst[1], β1_1_lst, 1, VWx_avg, VWy_avg, m, h, w, A, Cᴰ, Cᴸ, 0, 0,r_island];
package_2 = [α1_2_lst[1], α1_2_lst, β1_2_lst[1], β1_2_lst, 1, VWx_avg, VWy_avg, m, h, w, A, Cᴰ, Cᴸ, 0, 0,r_island];



# Define callbacks
cbs1, cbs2 = def_callbacks()

# Solve Equation of motion for boat 1

u1, t_1, t_β1, α1_1list, β1_1list = solve_boating(eq_motion, ut_o, u_o1,package_1, t_int,cbs1)


# Solve Equation of motion for boat 2
u2, t_2, t_β2, α1_2list, β1_2list = solve_boating(eq_motion, ut_o, u_o2,package_2, t_int,cbs2)


# Animate

anim_boat(u1,β1_1list,α1_1list,u2,β1_2list,α1_2list,VWx_avg,VWy_avg,VW_avg,U_avg,t_1,t_2,t_β1,t_β2,5)