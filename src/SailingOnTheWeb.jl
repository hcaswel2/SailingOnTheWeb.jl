module SailingOnTheWeb

using DifferentialEquations 
using LinearAlgebra  
using Interpolations
using ProgressMeter
using Plots
using Dierckx
using QuadGK
using Distributions, Random

export eq_motion, def_constants, generate_wind, def_callbacks, boatplot, anim_boat

function eq_motion(utt, ut, u, p, t)
    
    α1, α1_lst, β1, β1_lst, corner, VWx, VWy, m, h, w, Cᴰ, Cᴸ, corr_time, time, r_island = p 
    # I know this is alot... Many are just needed for callbacks
    
    V_boat = sqrt(ut[1]^2 + ut[2]^2)
    
    γ = atan(VWy(u[1],u[2],t[1])/VWx(u[1],u[2],t[1])) # Angle of aperent wind to the base wind velocity 
    # (island effects on wind angle)
    
    direct = sign(β1 + γ) # track if wind is hitting left side of boat or right side
    β1 = abs(β1) # Tread both sides the same using the variable "direct" to correct values
    
    VA = sqrt((ut[2]-VWy(u[1],u[2],t[1]))^2 + (VWx(u[1],u[2],t[1])-ut[1])^2) #m/s Apparent Wind speed
    θ = (γ + atan((ut[2]-VWy(u[1],u[2],t[1]))/(VWx(u[1],u[2],t[1])-ut[1])))*direct # Angle between wind and apparent wind
    β = direct*(β1 - γ) - θ # Angle of boat to apparent wind
    
    # Not all diagrams show this resistance done below, so I am not 100% sure this is not taken into account 
    # with lift and drag coefficients, but without some resistance, the boats 1 & 2 could catch Voyager 1 & 2
    # https://smalltridesign.com/Trimaran-Articles/Boat-Resistance.html
    Cs = 0.11
    Sp = (h * sin(α1)*w)/2
    Ra = 0.453592 * Cs * (VA*cos(β))^2 * Sp * -sign(ut[1])*5
    
    
    α = β - α1 # angle of sail to aparent wind
    ρ = 1.225 #kg/m^3 Density of air
    L = 1/2*ρ*VA^2*A*Cᴸ(α)
    D = 1/2*ρ*VA^2*A*Cᴰ(α)
    
    F_T = sqrt(L^2 + D^2) # find total force on boat
    F_LAT = -cos(β1-α1)*F_T #L*sin(α)-D*cos(α) #-cos(β1-α1)*F_T # find the amount pushing the boat to the side
    F_R = sin(β1-α1)*F_T #L*cos(α)-D*sin(α) #sin(β1-α1)*F_T # find the amount pushing the boat forward
    
    
    # Assume only 0% of the lateral force is contributing to the position
    utt[1,1] = (F_R-Ra)*-cos(β1)/m #+ 0.00*(F_LAT+Ra_LAT)*-sin(β1)/m #  # change to global coodinates x 
    utt[2,1] = (F_R-Ra)*sin(β1)/m*direct #- 0.00*(F_LAT+Ra_LAT)*cos(β1)/m*direct # change to global coodinates y
end

function def_constants()
    α_list = [10,20,30,40,50,60,70,80,90,100] .* pi ./ 180
    Cᴸ_list = [0.48,4.25,4.92,4.69,4.20,3.02,2.33,1.41,0.92,0.47] ./ 3.22
    Cᴰ_list = [0.48,0.58,0.71,0.95,1.30,1.78,2.42,3.14,4.10,5.15] ./ 3.22
    # https://en.wikipedia.org/wiki/Forces_on_sails

    w = (224/12/3.28084) #Sail width at bottom (assume triangle sail)
    h1 = 1.5 # Boom Height
    h2 = h1 + (8/3.28084) #top of sail height

    A = 1/2*(h2-h1)*w # Sail area
    h_avg = h1 + (h2-h1)/3 #average sail height

    r_island = 20; #island radius

  
    m = 59 + 65*2 #kg
    ρ = 1.225 #kg/m^3 air density

    

    # Create lift and drag coefficient arrays
    Cᴸ = Spline1D(α_list,Cᴸ_list)
    Cᴰ = Spline1D(α_list,Cᴰ_list)

    a = (10:1:100)
    CL = zeros(length(a))
    CD = zeros(length(a))
    for i in eachindex(a)
        CL[i] = Cᴸ(a[i]*pi/180)
        CD[i] = Cᴰ(a[i]*pi/180)
    end

    h = h2 - h1 # Sail height

    ut_o = [0.0; 0.0]
    t_int = 0.0:0.01:120

    return h, m, h, Cᴰ, Cᴸ, r_island, ut_o, t_int, A, h_avg, ρ
end

function generate_wind(h_avg,r_island)
    # Lets generate some random wind with a weibull distribution and a 
    # power spectral density function provided by Professor Shields

    kw = 2 # Shape param
    λw = 6 # Scale param
    X_cdf = 0:0.1:20
    cdf_W = cdf(Weibull(kw,λw),X_cdf) # Create the pdf for X2
    Weibull_invcdf = Spline1D(cdf_W,X_cdf)

    ω = 0:0.001:100
    t = 0:0.1:400 # Include negative time because the wind is measured from the center of the island and so 
    # upwind, the same base wind velocity is felt earlier 
    z = h1:(h2-h1)/10:h2

    # https://www.researchgate.net/publication/224327460_Peak_factor_estimation_in_hurricane_surface_winds
    Suu = zeros(length(z),length(ω)) 
    SP_funct = zeros(length(z),length(t))
    Weibull_wind = zeros(length(z),length(t))
    Weibull_mean_lst = zeros(length(z))

    dω = ω[2] - ω[1]

    for i in eachindex(z)
        Tu = 3.13 * z[i] ^0.2
        Suu[i,:] = 4 .* Tu ./ (1 .+ 70.8 .* (ω .* Tu).^2) .^(5/6) #Power specrtal density function shows wind 
        # relation over time
        
        # Perform the  spectral representation method to generate Random wind
        for j in eachindex(ω)
            V = rand(Normal(0, sqrt(2*Suu[i,j]*dω)), 1)
            W = rand(Normal(0, sqrt(2*Suu[i,j]*dω)), 1)
            SP_funct[i,:] += W.* cos.(ω[j].*t) + V.*sin.(ω[j].*t)
        end
        SP_funct[i,:] = SP_funct[i,:]
        Weibull_wind[i,:] = Weibull_invcdf(cdf(Normal(),SP_funct[i,:])) # Transfrom from Normal to Uniform and then to Wiebull
        Weibull_mean_lst[i] = mean(Weibull_wind[i,:])
    end
    mean_wind_funct = Spline1D(z,Weibull_mean_lst) # interpolate the wind at various height
    mean_wind = mean_wind_funct(h_avg) # Wind at the average sail height
    U_1 = Spline2D(z,t,Weibull_wind)
    U(z,t,x) = U_1(z,t-x/mean_wind)
    U_avg(t) = U(h_avg,t,0) 

    # Add island effects on wind
    VWx(z,x,y,t) = U(z,t,x) * (r_island^2*(y^2-x^2)+(x^2+y^2)^2)/(x^2+y^2)^2 #m/s Wind speed
    VWy(z,x,y,t) = -U(z,t,x) * (2*r_island^2*x*y)/(x^2+y^2)^2
    VW(z,x,y,t) = sqrt(VWx(z,x,y,t)^2 + VWy(z,x,y,t)^2)

    # windspeed at average sail height
    VWx_avg(x,y,t) = VWx(h_avg,x,y,t)
    VWy_avg(x,y,t) = VWy(h_avg,x,y,t)
    VW_avg(x,y,t) =  VW(h_avg,x,y,t);

    return VWx_avg, VWy_avg, VW_avg
end

function def_callbacks()
    function corner1(u,t,int) # Change direction (turn 1)
        int.u[3]+int.u[4]^2 < -int.p[15] && int.p[5] == 1
    end;
    
    function corner2(u,t,int) # Change direction (turn 2)
        -int.u[3]^2+int.u[4] > int.p[15] && int.p[5] == 2
    end;
    
    function corner3(u,t,int) # Change direction (turn 3)
        int.u[3]-int.u[4]^2 > int.p[15] && int.p[5] == 3
    end;
    
    function corner4(u,t,int) # Change direction (turn 4)
        int.u[3]^2+int.u[4] < -int.p[15] && int.p[5] == 4
    end;
    
    function enter_bounds1(u,t,int) # Correct to avoid island collision
        int.u[3]^2+int.u[4]^2 < (int.p[15]+8)^2 && int.p[14] > int.p[13]
        
    end;
    
    function exit_bounds1(u,t,int)  # Correct headding after correcting to avoid island collision
        int.u[3]^2+int.u[4]^2 >= (int.p[15]+8)^2 && int.p[14] > int.p[13] && int.p[3] != int.p[4][int.p[5]]
        
    end;
    
    function corner1_action!(int)
        display("Corner 1")
        angle = int.p[3]-int.p[4][2]
        int.p[1] = int.p[2][2]
        int.p[3] = int.p[4][2]
    
        if angle < pi/2 && angle > pi/2 # Depending on the turn angle, some velocity can be maintained
            x2 = cos(angle)*int.u[1] - sin(angle)*int.u[2]
            y2 = sin(angle)*int.u[1] + cos(angle)*int.u[2]
            int.u[1] = cos(angle)*x2
            int.u[2] = cos(angle)*y2
        else
            int.u[1] = 0
            int.u[2] = 0
        end
        int.p[5] = 2
    end;
    
    function corner2_action!(int)
        display("Corner 2")
        angle = int.p[3]-int.p[4][3]
        int.p[1] = int.p[2][3]
        int.p[3] = int.p[4][3]
        if angle < pi/2 && angle > pi/2 # Depending on the turn angle, some velocity can be maintained
            x2 = cos(angle)*int.u[1] - sin(angle)*int.u[2]
            y2 = sin(angle)*int.u[1] + cos(angle)*int.u[2]
            int.u[1] = cos(angle)*x2
            int.u[2] = cos(angle)*y2
        else
            int.u[1] = 0
            int.u[2] = 0
        end
        int.p[5] = 3
    end;
    
    function corner3_action!(int)
        display("Corner 3")
        angle = int.p[3]-int.p[4][4]
        int.p[1] = int.p[2][4]
        int.p[3] = int.p[4][4]
        
        if angle < pi/2 && angle > pi/2 # Depending on the turn angle, some velocity can be maintained
            x2 = cos(angle)*int.u[1] - sin(angle)*int.u[2]
            y2 = sin(angle)*int.u[1] + cos(angle)*int.u[2]
            int.u[1] = cos(angle)*x2
            int.u[2] = cos(angle)*y2
        else
            int.u[1] = 0
            int.u[2] = 0
        end
        int.p[5] = 4
    end;
    
    function corner4_action!(int)
        display("Corner 4")
        angle = int.p[3]-int.p[4][1]
        int.p[1] = int.p[2][1]
        int.p[3] = int.p[4][1]
        
        if angle < pi/2 && angle > pi/2 # Depending on the turn angle, some velocity can be maintained
            x2 = cos(angle)*int.u[1] - sin(angle)*int.u[2]
            y2 = sin(angle)*int.u[1] + cos(angle)*int.u[2]
            int.u[1] = (pi/2-angle)/pi*x2
            int.u[2] = (pi/2-angle)/pi*y2
        else
            int.u[1] = 0
            int.u[2] = 0
        end
        int.p[5] = 1
    end;
    
    function enter_bounds1_action!(int)
        display("Corrected course")
    
        int.p[3] = int.p[3] - 5 *pi/180
        angle = 5 *pi/180
        
        # Some velocity will be maintained
        x2 = cos(angle)*int.u[1] - sin(angle)*int.u[2]
        y2 = sin(angle)*int.u[1] + cos(angle)*int.u[2]
        int.u[1] = cos(angle)*x2
        int.u[2] = cos(angle)*y2
    
        int.p[13] = int.p[14] + 1 # allow time to pass before correcting again
    end;
    
    function exit_bounds1_action!(int)
        display("Recorrected course")
    
        int.p[3] = int.p[3] + 5 *pi/180
        angle = -5 *pi/180
        
        # Some velocity will be maintained
        x2 = cos(angle)*int.u[1] - sin(angle)*int.u[2]
        y2 = sin(angle)*int.u[1] + cos(angle)*int.u[2]
        int.u[1] = cos(angle)*x2
        int.u[2] = cos(angle)*y2
    
        int.p[13] = int.p[14] + 1 # allow time to pass before correcting again
    end;
    
    function saving_funct(u,t,int) # track the sail angle and boat angle
        int.p[14] = t 
        return [int.p[1],int.p[3]]
    end;

    # Define Callbacks
    cb1 = DiscreteCallback(corner1, corner1_action!);
    cb2 = DiscreteCallback(corner2, corner2_action!);
    cb3 = DiscreteCallback(corner3, corner3_action!);
    cb4 = DiscreteCallback(corner4, corner4_action!);
    cb_correction1 = DiscreteCallback(enter_bounds1, enter_bounds1_action!);
    cb_recorrection1 = DiscreteCallback(exit_bounds1, exit_bounds1_action!);

    saved_values1 = SavedValues(Float64, Vector{Float64})
    saving_cb1 = SavingCallback(saving_funct, saved_values1)

    saved_values2 = SavedValues(Float64, Vector{Float64})
    saving_cb2 = SavingCallback(saving_funct, saved_values2)

    cbs1 = CallbackSet(cb1, cb2, cb3, cb4, cb_correction1,cb_recorrection1, saving_cb1);
    cbs2 = CallbackSet(cb1, cb2, cb3, cb4, cb_correction1,cb_recorrection1, saving_cb2);
    return cbs1, cbs2
end

# This may be able to be rewritten more compactly, but I kept modifying so it got a little out of hand with 
# the input count
function boatplot(u,ux,uy,utx,uty,β1,α1,index,scale,i,γ,VW,l_color,plt_name)
    plot!([vcat(u[3][1:index],ux)],[vcat(u[4][1:index],uy)], linecolor = l_color, legendfontsize = 9, label = plt_name * string(round(sqrt(utx^2+uty^2),digits=1)) * " m/s")
    plot!([ux+cos(-β1)*scale,ux-cos(-β1)*scale],[uy+sin(-β1)*scale,uy-sin(-β1)*scale], linecolor = :black, label = false,linewidth=3) # Plot boat direction
    if β1 +γ < 0
        α1 = -α1
    end
    plot!([ux,ux+cos(-β1+α1)*scale*1.5],[uy,uy+sin(-β1+α1)*scale*1.5],linecolor = :black, label = false) # Plot sail direction
end

# Lets make an animator function
function anim_boat(u1,β1_1,α1_1,u2,β1_2,α1_2, VWx,VWy,VW, U, t_1, t_2,t_β1,t_β2,fps)
    # For real time we need a lot of spline/interpolation functions
    t_anim = (t_1[1]:1/fps:t_1[end])
    
    t1_index = LinearInterpolation(t_1,eachindex(t_1))
    t2_index = LinearInterpolation(t_2,eachindex(t_2))
    
    # Spline doesnt work because there are repeated values, so I used Interpolations.
    # Interpolations gets unhappy with the longest warning ever seen, but still runs, just need to clear the warning:
    IJulia.clear_output(true) 
    
    u1x = LinearInterpolation(t_1,u1[3])
    u1y = LinearInterpolation(t_1,u1[4])
    ut1x = LinearInterpolation(t_1,u1[1])
    ut1y = LinearInterpolation(t_1,u1[2])
    
    u2x = LinearInterpolation(t_2,u2[3])
    u2y = LinearInterpolation(t_2,u2[4])
    ut2x = LinearInterpolation(t_2,u2[1])
    ut2y = LinearInterpolation(t_2,u2[2])
    
    β11 = LinearInterpolation(t_β1,β1_1)
    α11 = LinearInterpolation(t_β1,α1_1)
    
    β12 = LinearInterpolation(t_β2,β1_2)
    α12 = LinearInterpolation(t_β2,α1_2)
    
    xmin = minimum(vcat(u1[3],u2[3]))
    xmax = maximum(vcat(u1[3],u2[3]))
    ymin = minimum(vcat(u1[4],u2[4]))
    ymax = maximum(vcat(u1[4],u2[4]))
    
    xlims = (xmin-(xmax-xmin)*0.1,xmax+(xmax-xmin)*0.1)
    ylims = (ymin-(ymax-ymin)*0.1,ymax+(ymax-ymin)*0.1 + 10)
    scale = (max(xmax,ymax) - min(xmin,ymin))/40

    t_anim_end = eachindex(t_anim)[end]
    prog = Progress(t_anim_end, "Building Animation: ")
    
    anim = @animate for i in eachindex(t_anim)
        plot(size = (800,800),xlim =xlims,ylim=ylims,aspect_ratio=:equal,legend = :topright,framestyle=:box,title = "Base Wind Velocity = " * string(round(U(t_anim[i]),digits=1)))
        index1 = Int(floor(t1_index(t_anim[i])))
        index2 = Int(floor(t2_index(t_anim[i])))
        xs = xlims[1]:(xlims[end]-xlims[1])/20:xlims[end]
        ys = ylims[1]:(ylims[end]-ylims[1])/20:ylims[end]

        df(x, y) = normalize([VWx(x,y,t_anim[i]), VWy(x,y,t_anim[i])]) .*3

        xxs = [x for x in xs for y in ys if x^2+y^2 > r_island^2] # Dont Plot in the island
        yys = [y for x in xs for y in ys if x^2+y^2 > r_island^2] # Dont Plot in the island

        quiver!(xxs, yys, quiver=df) #Vector field generation
        
        k = 0:pi/100:2*pi
        x1a = sin.(k)*(r_island+8)
        y1a = cos.(k)*(r_island+8)
        x2 = sin.(k)*r_island
        y2 = cos.(k)*r_island
        
        γ1 = atan(VWy(u1x(t_anim[i]),u1y(t_anim[i]),t_anim[i])/VWx(u1x(t_anim[i]),u1y(t_anim[i]),t_anim[i]))
        γ2 = atan(VWy(u2x(t_anim[i]),u2y(t_anim[i]),t_anim[i])/VWx(u2x(t_anim[i]),u2y(t_anim[i]),t_anim[i]))
        
        plot!(x1a,y1a,linecolor = :black,linestyle = :dash, label = false)
        plot!(x2,y2, label = false,seriestype = [:shape,],fillcolor = :sienna4)
        boatplot(u1,u1x(t_anim[i]),u1y(t_anim[i]),ut1x(t_anim[i]),ut1y(t_anim[i]),β11(t_anim[i]),α11(t_anim[i]),
            index1,scale,i, γ1,VW(u1x(t_anim[i]),u1y(t_anim[i]),t_anim[i]),:red, "Boat 1 Velocity = ")
        boatplot(u2,u2x(t_anim[i]),u2y(t_anim[i]),ut2x(t_anim[i]),ut2y(t_anim[i]),β12(t_anim[i]),α12(t_anim[i]),
            index2,scale,i, γ2,VW(u2x(t_anim[i]),u2y(t_anim[i]),t_anim[i]),:blue, "Boat 2 Velocity = ")

        next!(prog)
    end
    gif(anim, "Boat1.mp4",fps = fps)
end;

end # module



