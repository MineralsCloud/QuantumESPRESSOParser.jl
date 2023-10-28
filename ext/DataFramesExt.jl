module DataFramesExt

using DataFrames: DataFrame, vcat
using QuantumESPRESSOParser.PWscf: eachstep, eachconvergedenergy, eachunconvergedenergy

import QuantumESPRESSOParser.PWscf: each_energy_by_step

function each_energy_by_step(str::AbstractString)
    return Iterators.map(
        zip(eachstep(str), eachconvergedenergy(str))
    ) do (str_at_step, converged)
        unconverged_df = DataFrame(collect(eachunconvergedenergy(str_at_step)))
        converged_df = DataFrame(Base.vect(converged))  # It's only one row.
        vcat(unconverged_df, converged_df; cols=:union)
    end
end

end
