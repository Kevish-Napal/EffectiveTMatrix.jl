## ParticulateSphere is a specific alias of the type Material defined in EffectiveWaves.jl

ParticulateSphere{T,Dim,P,Sps} = Material{Sphere{T,Dim},M} where M <: ParticulateMicrostructure{Dim,P,Sps}
ParticulateCylinder{T} = ParticulateSphere{T,2}

radius(sphere::ParticulateSphere) = sphere.shape.radius;



# struct EffectiveParticle{T,M<:Microstructure} where T
#     radius::T
#     microstructure::M
#     function EffectiveParticle{M}(radius::T,micro::M) where {T,M <: Microstructure}
#         new{T,M}(radius,micro)
#     end
# end



# import EffectiveWaves: shape

# function shape(p::EffectiveParticle{T,Dim}) where {T,Dim}
#     Sphere{T,Dim}(p.radius) 
# end



function t_matrix(ω::AbstractFloat, host_medium::PhysicalMedium, material::ParticulateSphere; basis_order=10::Integer, basis_field_order=10::Integer)

    kstar, wavemode = solve_eigensystem(ω, host_medium, material; basis_order=basis_order);
    N,D = t_matrix_num_denom(kstar,wavemode;basis_field_order=basis_field_order);
    T =  (- vec(sum(N,dims=1)./sum(D,dims=1)))
    # T0 = (- N[1+basis_order,:]./D[1+basis_order,:])

    return T #, T0
end


function solve_eigensystem(ω::AbstractFloat, 
    host_medium::Acoustic{T,dim}, material::ParticulateSphere{T,dim}; basis_order=10::Integer) where {T,dim}
    
    sp_EF = sp_MC_to_EF(material.microstructure.species[1],material.shape.radius)
    
    opts = Dict(
       :tol => 1e-4, 
       :num_wavenumbers => 1
       ,:basis_order => basis_order
   );
   
   kstar = wavenumbers(ω,host_medium, sp_EF;opts...)[1]

   kws = Dict(
       :basis_order => basis_order
       ,:tol=>1e-2
   );
   
   rsource = regular_spherical_source(host_medium, [1.0+0.0im];
   position = zeros(dim), symmetry = RadialSymmetry{dim}()
   );

   wavemode = WaveMode(ω, kstar, rsource, material; kws...);

   return kstar, wavemode
end 

function t_matrix_num_denom(kstar,wavemode;basis_field_order)
    
    material = wavemode.material
    micro = material.microstructure
    
    k=wavemode.ω/micro.medium.c
    R = material.shape.radius
    species = micro.species
    Rtildas = R .- outer_radius.(species)

    F = wavemode.eigenvectors
    nbo, n_λ, _ = size(F)
    basis_order = Int((nbo-1)/2)
    L = basis_order+basis_field_order


    n_densities = [number_density(sp) for sp in species]

    J = [besselj(n,k*Rtildas[λ])*n_densities[λ] 
        for n = -L-1:L,
            λ in 1:n_λ]

    Jstar = [besselj(n,kstar*Rtildas[λ])*n_densities[λ] 
        for n = -L-1:L,
            λ in 1:n_λ]

    H = [besselh(n,k*Rtildas[λ])*n_densities[λ] 
        for n = -L-1:L,
            λ in 1:n_λ]

    pre_num = k*J[1:1+2*L,:].*Jstar[2:2*(1+L),:] -
              kstar*J[2:2*(1+L),:].*Jstar[1:1+2*L,:]
    pre_denom = k*H[1:1+2*L,:].*Jstar[2:2*(1+L),:] -
                kstar*H[2:2*(1+L),:].*Jstar[1:1+2*L,:]

    Matrix_Num = complex(zeros(2*basis_order+1,2*basis_field_order+1))
    Matrix_Denom = complex(zeros(2*basis_order+1,2*basis_field_order+1))

    for n = 1:1+2*basis_field_order
        
        Matrix_Num[:,n] = vec(sum(pre_num[n+2*basis_order:-1:n,:].*F,dims=2))

        Matrix_Denom[:,n] = vec(sum(pre_denom[n+2*basis_order:-1:n,:].*F,dims=2))
    end

    return Matrix_Num, Matrix_Denom
end

function average_scattered_field(ω::AbstractFloat, region, source::AbstractSource{P}, material::ParticulateSphere{T,Dim};
    res=300, kws...) where {T,Dim,P<:PhysicalMedium{Dim},V <: AbstractArray{T}}

    x_vec, _ = points_in_shape(region;resolution=res);
    average_scattered_field(ω, x_vec, source, material; kws...)
end

function average_scattered_field(ω::AbstractFloat, x_vec::Vector{V}, source::AbstractSource{P}, material::ParticulateSphere{T,Dim};
    only_scattered_waves=false::bool, basis_order=10::Integer, basis_field_order=10::Integer) where {T,Dim,P<:PhysicalMedium{Dim},V <: AbstractArray{T}}

    T_matrix =  t_matrix(ω, source.medium, material; basis_order=basis_order, basis_field_order=basis_field_order)
    source_coefficient = regular_spherical_coefficients(source)
    G = source_coefficient(basis_field_order,zeros(Dim),ω)
    
    N = basislength_to_basisorder(P,length(G))
    basis = outgoing_basis_function(source.medium, ω)
    
    function us_func(x)
        basis(N, x)*(G.* T_matrix)
    end
    
    # println("T_matrix = ", T_matrix)
    # println("G = ", G)
    # print("basis=",basis(N, [.324,0.0]))
    # print("combination=",basis(N, [.324,0.0])*(G.* T_matrix))

    if only_scattered_waves
        us = map(x_vec) do x
            if sqrt(x'*x) > radius(material)
                us_func(x)
            else
                .0+.0im # unknown for now
            end
        end
    
        return FrequencySimulationResult(us,x_vec,[ω])

    else
        utot = map(x_vec) do x
            if sqrt(x'*x) > radius(material)
                field(source)(x,ω) + us_func(x)
            else
                .0+.0im # unknown for now
            end
        end
    
        return FrequencySimulationResult(utot,x_vec,[ω])
    end
end