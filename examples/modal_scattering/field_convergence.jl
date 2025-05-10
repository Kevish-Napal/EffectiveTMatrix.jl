
using Pkg
Pkg.activate("EffectiveTMatrix")
using EffectiveTMatrix
using Plots
using Statistics

dimension=2;                                       
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);   # 2D acoustic problem
input_mode = 3;                                          # mode n=0 
source =   mode_source(host_medium,input_mode); # V_0


ω=0.1; #0.2                                    # frequency
M=N=45.0;  
resolution = 200                                 # sizes of the rectangle where to plot
bottomleft = [-M;-N]; topright = [M;N];
region = Box([bottomleft, topright]);



# parameters for particles configurations
radius_big_cylinder = 20.0;                                 # radius of cylinder where particles are confined 
particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0)); # particles type: sound hard particles of radius 1
sp_MC = Specie(particle; volume_fraction = .05);            # # volume fraction is the density of particles

microstructure = Microstructure(host_medium,[sp_MC_to_EF(sp_MC, radius_big_cylinder)]); # microstructure of the medium
cylinder = Material(Circle(radius_big_cylinder),microstructure);
us = average_scattered_field(ω, region, source, cylinder; only_scattered_waves=true,basis_field_order=5);
# plot(us, ω; field_apply=real,seriestype = :contour,c=:balance)
expected_plot = plot(us, ω; field_apply=real,seriestype = :contour,c=:balance,legend=:none,title="Predicted\nEnsemble Average\n");
expected_plot = annotate!(0.0,0.0,text("effective\nmaterial", :black, 10))

basis_order=5;

x_vec, _ = points_in_shape(region;resolution=resolution);                            # space discretization 
nb_of_configurations = 100;                                                   # total number of config 

l = @layout [a b c; d{0.1h}]

A = complex(zeros(length(x_vec),nb_of_configurations));                       # store the fields of each config
anim = @animate for i ∈ 1:nb_of_configurations

    # plot scattering field
    particles_realisation = renew_particle_configurations(sp_MC,radius_big_cylinder);
    sim = FrequencySimulation(particles_realisation,source);
    scattered_field = run(sim,region,[ω];only_scattered_waves=true,basis_order=basis_order,res=resolution);
    sc_plot = plot(scattered_field,ω; field_apply=real,seriestype = :contour,c=:balance,title="\n Single Realisation\n",legend=:none);
    
    # plot average scattered field
    A[:,i] = mean.(scattered_field.field[:,1])
    mean_A  = mean(A[:,1:i],dims=2);
    mean_us = FrequencySimulationResult(mean_A,x_vec,[ω]);
    mean_sc_plot = plot(mean_us,ω; field_apply=real,seriestype = :contour,c=:balance,title="Cumulative\nEnsemble Average\n",legend=:none);

    
    # Progress bar dimensions (bar background and fill)
    progress = i / nb_of_configurations
    bar_y = 0.0
    progress_plot = plot([0, 2π], [bar_y, bar_y], lw=6, color=:lightgray, label = "", grid=false,axis=false, ylims=(-.01,.01),size=(1500, 100))
    progress_plot = plot!([0, 2π * progress], [bar_y, bar_y], lw=6, color=:gray, label = "")
    
    plot(sc_plot,mean_sc_plot, expected_plot, progress_plot,layout=l,size=(1500,700))
end
gif(anim, "animation1.gif", fps = nb_of_configurations/100)



