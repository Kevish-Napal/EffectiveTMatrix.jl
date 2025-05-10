# Return matrix of computed coefficients a col corresponds to a configuration

function mode_source(medium::Acoustic{T,2},N::Int) where T <: AbstractFloat
    coeffs = complex(zeros(2abs(N)+1))
    N<0 ? coeffs[1] = 1.0+0.0im : coeffs[end] = 1.0+0.0im
    return regular_spherical_source(medium, coeffs; 
        position = zeros(2),
        symmetry = (N==0) ? RadialSymmetry{2}() : WithoutSymmetry{2}()
        );
end

# this function handles Spheres containing only one type of specie 
function renew_particle_configurations(sp::Specie,radius_big_cylinder::Float64)
    config = random_particles(
        sp.particle.medium,
        sp.particle.shape;
        region_shape = Circle(1.05*radius_big_cylinder),
        volume_fraction = sp.volume_fraction
    );
    config = config[norm.(origin.(config)) .< radius_big_cylinder .- outer_radius.(config)]
end


function run_MC_validation!(host_medium::PhysicalMedium{Dim}, MC::MonteCarloResult{N};kws...) where {Dim,N}
    run_MC_validation!(host_medium, [MC]; kws...);
    return nothing
end

# main function for run_MC_validation
function run_MC_validation!(host_medium::PhysicalMedium{Dim}, MC_vec::Vector{MonteCarloResult{N}};
        basis_order::Int=10, kws...) where {Dim,N}

    basis_field_order = N-1
    for MC in MC_vec      
        # Effective method
        micro = Microstructure(host_medium,[MC.sp_MC]);
        material = Material(Sphere{Float64,Dim}(zeros(Dim),MC.R),micro);
        μeff = t_matrix(MC.ω, host_medium, material, basis_order=basis_order, basis_field_order=basis_field_order)[basis_field_order+1:2basis_field_order+1];
        # Monopole approximation of the T matrix
        μeff0 = t_matrix(MC.ω, host_medium, material, basis_order=0, basis_field_order=basis_field_order)[basis_field_order+1:2basis_field_order+1];
        MC.μeff = μeff
        MC.μeff0 = μeff0
    end

    return nothing
end


function update_monopole_approx!(host_medium::PhysicalMedium{Dim}, MC_vec::Vector{MonteCarloResult{N}}) where {Dim,N}
    basis_field_order = N-1
    for MC in MC_vec
        micro = Microstructure(host_medium,[MC.sp_MC]);
        material = Material(Sphere{Float64,Dim}(zeros(Dim),MC.R),micro);
        μeff0 = t_matrix(MC.ω, host_medium, material, basis_order=0, basis_field_order=basis_field_order)[basis_field_order+1:2basis_field_order+1];
        MC.μeff0 = μeff0
    end
    return nothing
end

"""
    sample_effective_t_matrix(ω::Number, host_medium::PhysicalMedium, sp::Specie, select_modes::Vector{Bool};kws...)

Monte Carlo computation of the effective cylinder T_matrix. The coefficients are computed for n 
in [0; basis_field_order]
"""

function sample_effective_t_matrix!(MC_vec::Vector{MonteCarloResult{N}},host_medium::PhysicalMedium{Dim}; 
    kws...) where {N,Dim}

    L = length(MC_vec)
    l = 0
    for MC in MC_vec
        l+=1
        println(l,"/",L)
        sample_effective_t_matrix!(MC, host_medium; kws...)
    end
end

# This function is basically a for loop on parameter input_mode appearing in function optimal2_mode_analysis 
# with elementary additional optimization. For each configuration of particles,
# all the required bessel functions are computed at once in blocs_V. The matrix
# V appearing in the function optimal2_mode_analysis selects the appropriate bloc of blocs_V according to the mode computed.
# Finally the stop criteria of computing each modes are set independently (this is why basis_field_order is updated after 
# each loops). I still need to check if this is required or not, for the moment it is difficult to set a convergence criteria. 
# note that a mode is assumed to be zero as soon as smaller than 1e-8, the code then stops iterating on this mode.
# Monte Carlo Simulation for the acoustic cylinder
function sample_effective_t_matrix!(MC::MonteCarloResult{N},host_medium::PhysicalMedium{Dim};
    basis_order::Int, nb_iterations_max=500::Int,nb_iterations_step=200::Int,prec=1e-1::Float64) where {N,Dim}

    ω = MC.ω
    sps_MC = [MC.sp_MC]
    radius_big_cylinder = MC.R
    basis_field_order = N-1

    
    k = ω/host_medium.c
    # Precompute T-matrices for these particles
    TMatrix = get_t_matrices(host_medium, [sp.particle for sp in sps_MC], ω, basis_order)[1]

    # We need this big vector F to store realisation of each modes 
    # which convergent rate my differ
    F = [ComplexF64[] for _ in 0:basis_field_order]

    # To estimate convergence of each modes, we keep track of number of iterations of each of them:
    total_iterations = 0;  

    initial_basis_field_order = basis_field_order
    select_modes = trues(basis_field_order+1)

    
    
    # confidence interval level 95% is [m-1.96σ/√n ; m+1.96σ/√n]
    # we have [1.96σ/√nb_iterations < prec*|m|] => m = empirical_mean ± prec*|m|] 
    # mode_continue_crit(m,s_r,s_i,n,prec) = (1.96*s_r/sqrt(n) > abs(real(m))*prec || 1.96*s_i/sqrt(n) > abs(imag(m))*prec) && abs(m) > 1e-8
    # the function below updates select_modes accordingly
    # function update_running_modes!(F,basis_field_order,select_modes,total_iterations,prec)
    #     for mode=0:basis_field_order
    #         if select_modes[mode+1]
    #             m = mean(F[mode+1]) # can be computed iteratively and stored to optimize
    #             s_r = std(real.(F[mode+1]); mean=real(m))
    #             s_i = std(imag.(F[mode+1]); mean=imag(m))
    #             select_modes[mode+1] = (1.96*s_r/sqrt(total_iterations) > abs(real(m))*prec || 1.96*s_i/sqrt(total_iterations) > abs(imag(m))*prec) && abs(m) > 1e-8
    #         end
    #     end
    #     if any(select_modes)
    #         basis_field_order = maximum(collect(0:initial_basis_field_order)[select_modes])
    #     end
    # end

    while any(select_modes) && maximum(total_iterations) < nb_iterations_max
        for _ in 1:nb_iterations_step

            # note sps_MC[1] in renew_particle_configurations - the latter function only takes into account one type of sp for now
            particles = renew_particle_configurations(sps_MC[1],radius_big_cylinder)
            rθ = [cartesian_to_radial_coordinates(origin(p)) for p in particles]
            n_particles = length(particles)

            # M_bessel = [Jₙ(krᵢ)exp(Im*n*θᵢ) for n=0:bo+bfo, i=1:nb_particles] (no mistake: n starts at 0)
            blocs_V₊ = [ besselj(n,k*rθ[i][1])*exp(im*n*rθ[i][2]) 
                                    for n=0:basis_order+basis_field_order, i=1:n_particles]   
            
            blocs_V₋ = [ (-1)^n*conj(blocs_V₊[n+1,i])
                                    for n=1:basis_order, i=1:n_particles] 

            blocs_V = reverse(vcat(reverse(blocs_V₋,dims=1),blocs_V₊),dims=1)
            # blocs_V is of size (2*basis_order + basis_field_order + 1) x n_particles

            # Compute scattering matrix for all particles
            S = scattering_matrix(host_medium, particles, [TMatrix for p in particles], ω, basis_order)

            V = Array{ComplexF64,2}(undef,2*basis_order+1,n_particles)
            a = Array{ComplexF64,2}(undef,2*basis_order+1,n_particles)

      
            for input_mode = 0:basis_field_order # optimal when this loop comes after renewing particles
                input_mode_index = input_mode+1
                if  select_modes[input_mode_index]

                    inds = (basis_field_order-input_mode+1):(basis_field_order-input_mode+2*basis_order+1)
                    V .= blocs_V[inds,:]                  
                    
                    # reshape and multiply by t-matrix to get the scattering coefficients
                    a .= reshape((S + I) \ reduce(vcat,V),2*basis_order+1,n_particles)
                    for i in axes(a,2)
                        a[:,i] = TMatrix * a[:,i]
                    end

                    # this uses a lot of memory, should be optimized
                    push!(F[input_mode_index],sum(conj(V).*a))
                    
                end
            end # mode loop                                                      
        end # iteration step
       
        total_iterations += nb_iterations_step

        #update_running_modes!(F,basis_field_order,select_modes,total_iterations,prec)
        for mode=0:basis_field_order
            if select_modes[mode+1]
                m = mean(F[mode+1]) # can be computed iteratively and stored to optimize
                s_r = std(real.(F[mode+1]); mean=real(m))
                s_i = std(imag.(F[mode+1]); mean=imag(m))
                select_modes[mode+1] = (1.96*s_r/sqrt(total_iterations) > abs(real(m))*prec || 1.96*s_i/sqrt(total_iterations) > abs(imag(m))*prec) && abs(m) > 1e-8
            end
        end
        if any(select_modes)
            basis_field_order = maximum(collect(0:initial_basis_field_order)[select_modes])
        end
        ## finished updating running modes

        println("nb iterations:",total_iterations)
        println("modes still running:", select_modes,"\n")
    end # while
    
    MC.μ = mean.(F);
    MC.σ = [std(real.(F[N]);mean=real(MC.μ[N])) + im*std(imag.(F[N]);mean=imag(MC.μ[N])) for N=1:initial_basis_field_order+1];
    MC.nb_iterations = length.(F);
    return nothing
end


function naive_sample_effective_t_matrix(ω::Number, host_medium::PhysicalMedium, sp::Specie;
    radius_big_cylinder=10.0::Float64, basis_order=3::Int, basis_field_order=0::Int,
    nb_iterations=1000::Int) 

    k = ω/host_medium.c
    F = [ComplexF64[] for _ in 1:basis_field_order+1]
    x = 1.5*radius_big_cylinder
   
    for mode=0:basis_field_order
       source = mode_source(host_medium,mode)
       for _ = 1:nb_iterations
        particles = renew_particle_configurations(sp,radius_big_cylinder)
        sim = FrequencySimulation(particles,source);
        result = run(sim,[[x,0.0]],[ω];only_scattered_waves=true,basis_order=basis_order)
        Fnn = result.field[1][1]/besselh(mode,k*x)
        push!(F[mode+1],Fnn)
       end
    end 
    return F
end

# In this function we analyse the mode of the scattered field for given mode source order N
# outputs all realisation, not the average.
# This function is not optimized, it just checks the theory: all modes of <us> different from N should vanish
function mode_analysis(input_mode::Int, ω::Number, host_medium::PhysicalMedium, sp::Specie;
    radius_big_cylinder=10.0::Float64, basis_order=3::Int, basis_field_order=0::Int,
    nb_iterations=1000::Int) 

    k = ω/host_medium.c
    source = mode_source(host_medium,input_mode)

    function kernel_V(n,particles)
        V = complex(zeros(2*basis_order+1,length(particles)))
        for j = 1:length(particles)
            r, θ  = cartesian_to_radial_coordinates(-origin(particles[j]))
            V[:,j] = [besselj(m,k*r)*exp(im*θ*m) for m = -basis_order-n:basis_order-n]
        end
        return V
    end

    F0 = Array{ComplexF64}(undef,basis_field_order+1,nb_iterations)
    for iter = 1:nb_iterations
        particles = renew_particle_configurations(sp,radius_big_cylinder)
        sim = FrequencySimulation(particles,source);
        scattering_coefficients = basis_coefficients(sim, ω; basis_order = basis_order)
        F0[:,iter] = [sum(kernel_V(n,particles).*scattering_coefficients) for n = 0:basis_field_order]
    end
    
    return [F0[n,:] for n = 1:basis_field_order+1]
end

# Very similar to the function above, however we only compute the mode that is not supposed to vanish.
# Also, we prepare for the computation of V in a more optimal way
# (allowing computation of bessel functions least possible times)
function optimal1_mode_analysis(input_mode::Int, ω::Number, host_medium::PhysicalMedium, sp::Specie;
    radius_big_cylinder=10.0::Float64, basis_order=3::Int, nb_iterations=1000::Int) 

    k = ω/host_medium.c
    source = mode_source(host_medium,input_mode)

    F = ComplexF64[]
    for _ = 1:nb_iterations
        particles = renew_particle_configurations(sp,radius_big_cylinder)
        rθ = [cartesian_to_radial_coordinates(origin(p)) for p in particles]
        n_particles = length(particles)
        V = [besselj(n,k*rθ[i][1])*exp(im*n*rθ[i][2]) 
                for n = (input_mode+basis_order):-1:(input_mode-basis_order), i=1:n_particles]

        sim = FrequencySimulation(particles,source);
        scattering_coefficients = basis_coefficients(sim, ω; basis_order = basis_order)
        push!(F,sum(conj(V).*scattering_coefficients))
    end
    
    return F
end

# Very similar to the function above, however we optimize the simulation of the scattering problem
# For example functions mode_source, scattering_coefficients are not called (they would all recompute bessel functions)
function optimal2_mode_analysis(input_mode::Int, ω::Number, host_medium::PhysicalMedium, sp::Specie;
    radius_big_cylinder=10.0::Float64, basis_order=3::Int, basis_field_order=0::Int,
    nb_iterations=1000::Int) 

    t_matrix = get_t_matrices(host_medium, [sp.particle], ω, basis_order)[1]
    k = ω/host_medium.c
    # source = mode_source(input_mode)

    F = ComplexF64[]
    for _ = 1:nb_iterations
        particles = renew_particle_configurations(sp,radius_big_cylinder)
        rθ = [cartesian_to_radial_coordinates(origin(p)) for p in particles]
        n_particles = length(particles)
        V = [besselj(n,k*rθ[i][1])*exp(im*n*rθ[i][2]) 
                for n = (input_mode+basis_order):-1:(input_mode-basis_order), i=1:n_particles]

        
        # Compute scattering matrix for all particles
        S = scattering_matrix(host_medium, particles, [t_matrix for p in particles], ω, basis_order)

        # reshape and multiply by t-matrix to get the scattering coefficients
        a = reshape((S + I) \ reduce(vcat,V),2*basis_order+1,n_particles)
        for i in axes(a,2)
            a[:,i] = t_matrix * a[:,i]
        end
        push!(F,sum(conj(V).*a))
    end
    
    return F
end

function sp_MC_to_EF(sp_MC::Specie,radius_big_cylinder::Float64)

    mean_nb_particles = mean([length(renew_particle_configurations(sp_MC,radius_big_cylinder)) for i = 1:500])
    particle_radius = outer_radius(sp_MC.particle)
    ϕ = mean_nb_particles*particle_radius^2/(radius_big_cylinder-particle_radius)^2
    
    return Specie(Acoustic(2; ρ=sp_MC.particle.medium.ρ, c=sp_MC.particle.medium.c),Circle(particle_radius);
    volume_fraction = ϕ,separation_ratio=sp_MC.separation_ratio)
end

function generate_species(radius_big_cylinder::Float64,particle::Particle,ϕ::Float64,separation_ratio::Float64)

    sp_MC = Specie(particle; volume_fraction = ϕ,separation_ratio=separation_ratio) 
    sp_EF = sp_MC_to_EF(sp_MC,radius_big_cylinder)

    return sp_MC, sp_EF
end


function sample_effective_t_matrix_archive(ω::Number, host_medium::PhysicalMedium, sp::Specie;
    radius_big_cylinder=10.0::Float64, basis_order=10::Int, basis_field_order=0::Int,
    nb_iterations_max=5000::Int,nb_iterations_step=200::Int,prec=1e-1::Float64) 

    k = ω/host_medium.c
    # Precompute T-matrices for these particles
    t_matrix = get_t_matrices(host_medium, [sp.particle], ω, basis_order)[1]

    # We need this big vector F to store realisation of each modes 
    # which convergent rate my differ
    F = [ComplexF64[] for _ in 0:basis_field_order]

    # To estimate convergence of each modes, we keep track of number of iterations of each of them:
    total_iterations = 0;  

    initial_basis_field_order = basis_field_order
    select_modes = trues(basis_field_order+1)

    continue_crit() = any(select_modes) && maximum(total_iterations) < nb_iterations_max
    
    # confidence interval level 95% is [m-1.96σ/√n ; m+1.96σ/√n]
    # we have [1.96σ/√nb_iterations < prec*|m|] => m = empirical_mean ± prec*|m|] 
    # mode_continue_crit(m,s_r,s_i,n,prec) = (1.96*s_r/sqrt(n) > abs(real(m))*prec || 1.96*s_i/sqrt(n) > abs(imag(m))*prec) && abs(m) > 1e-8
    # the function below updates select_modes accordingly
    function update_running_modes!()
        for mode=0:basis_field_order
            if select_modes[mode+1]
                m = mean(F[mode+1]) # can be computed iteratively and stored to optimize
                s_r = std(real.(F[mode+1]); mean=real(m))
                s_i = std(imag.(F[mode+1]); mean=imag(m))
                select_modes[mode+1] = (1.96*s_r/sqrt(total_iterations) > abs(real(m))*prec || 1.96*s_i/sqrt(total_iterations) > abs(imag(m))*prec) && abs(m) > 1e-8
            end
        end
        if any(select_modes)
            basis_field_order = maximum(collect(0:initial_basis_field_order)[select_modes])
        end
    end

    while continue_crit()
        for _ in 1:nb_iterations_step

            particles = renew_particle_configurations(sp,radius_big_cylinder)
            rθ = [cartesian_to_radial_coordinates(origin(p)) for p in particles]
            n_particles = length(particles)

            # J = [Jₙ(krᵢ) for n=0:bo+bfo, i=1:nb_particles]
            J = [
                    besselj(n,k*rθ[i][1]) 
                        for n = 0:(basis_order+basis_field_order),i=1:n_particles
                ]

            # M_bessel = [Jₙ(krᵢ)exp(Im*n*θᵢ) for n=0:bo+bfo, i=1:nb_particles] (no mistake: n starts at 0)
            blocs_V₊ = [ besselj(n,k*rθ[i][1])*exp(im*n*rθ[i][2]) 
                                    for n=0:basis_order+basis_field_order, i=1:n_particles]   
            
            blocs_V₋ = [ (-1)^n*conj(blocs_V₊[n+1,i])
                                    for n=1:basis_order, i=1:n_particles] 

            blocs_V = reverse(vcat(reverse(blocs_V₋,dims=1),blocs_V₊),dims=1)
            # blocs_V is of size (2*basis_order + basis_field_order + 1) x n_particles

            # Compute scattering matrix for all particles
            S = scattering_matrix(host_medium, particles, [t_matrix for p in particles], ω, basis_order)

            V = Array{ComplexF64,2}(undef,2*basis_order+1,n_particles)
            a = Array{ComplexF64,2}(undef,2*basis_order+1,n_particles)

      
            for input_mode = 0:basis_field_order # optimal when this loop comes after renewing particles
                input_mode_index = input_mode+1
                if  select_modes[input_mode_index]

                    inds = (basis_field_order-input_mode+1):(basis_field_order-input_mode+2*basis_order+1)
                    V .= blocs_V[inds,:]                  
                    
                    # reshape and multiply by t-matrix to get the scattering coefficients
                    a .= reshape((S + I) \ reduce(vcat,V),2*basis_order+1,n_particles)
                    for i in axes(a,2)
                        a[:,i] = t_matrix * a[:,i]
                    end

                    # this uses a lot of memory, should be optimized
                    push!(F[input_mode_index], sum(conj(V).*a))
                    
                end
            end # mode loop                                                      
        end # iteration step
       
        total_iterations += nb_iterations_step
        update_running_modes!()

        println("nb iterations:",total_iterations)
        println("modes still running:", select_modes,"\n")
    end # while

    return F
end