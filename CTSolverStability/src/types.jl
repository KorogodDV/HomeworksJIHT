abstract type Component end

abstract type Mixture end

ncomponents(m::Mixture) == error("Not Implemented")

log_fugacity_coef(m::Mixture, molfrac, pressure, RT, phase::Symbol) = error("Not Implemented")