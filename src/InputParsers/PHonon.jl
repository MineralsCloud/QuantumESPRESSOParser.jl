module PHonon

using MLStyle: @match

using QuantumESPRESSOBase.Namelists: Namelist, to_dict
using Compat: isnothing
using Fortran90Namelists.FortranToJulia: FortranData
using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PHonon
using QuantumESPRESSOBase
using QuantumESPRESSOBase.Cards
using QuantumESPRESSOBase.Cards.PHonon: SpecialQPoint, 
                                        QPointsSpecsCard
using QuantumESPRESSOBase.Inputs
using QuantumESPRESSOBase.Inputs.PHonon

using QuantumESPRESSOParsers.InputParsers.Namelists

const Q_POINTS_SPECIAL_BLOCK_REGEX = r"""
^ [ \t]* qPointsSpecs [ \t]*$\n
^ [ \t]* \S+ [ \t]* $\n  # nqs
(?P<block>
 (?:
  ^ [ \t]* \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]+ \S+ [ \t]* $\n?
 )+
)
"""imx

const Q_POINTS_SPECIAL_ITEM_REGEX = r"""
^ [ \t]* (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]+ (\S+) [ \t]* $\n?
"""mx

function Base.parse(::Type{<:QPointsSpecsCard}, str::AbstractString)
    m = match(Q_POINTS_SPECIAL_BLOCK_REGEX, str)
    if !isnothing(m)
        captured = m.captures[1]
        data = SpecialQPoint[]
        for matched in eachmatch(Q_POINTS_SPECIAL_ITEM_REGEX, captured)
            # TODO: Match `nqs`
            point = @match map(x -> parse(Float64, FortranData(x)), matched.captures) begin
                [coordinates..., weight] => SpecialQPoint(coordinates, weight)
            end
            push!(data, point)
        end
        return QPointsSpecsCard(data)
    end

    @info "Cannot find card `Q_POINTS`!"
end # function Base.parse

function Base.parse(::Type{PHononInput}, str::AbstractString)
    dict = Dict{Symbol,InputEntry}()
    dict[name(PHNamelist)] = parse(PHNamelist, str)
    dict[:q_points] = parse(QPointsSpecsCard, str)
    return PHononInput(dict...)
end # function Base.parse

content=raw"""


&inputph

outdir='./tmp'
prefix='MgSiO3'
epsil=.true.
fildyn='dynmat'
tr2_ph=1.0d-14
ldisp=.true.
amass(1)=24.30500
amass(2)=28.08600
amass(3)=15.99900
nq1=4,nq2=4,nq3=4
/
qPointsSpecs
3
1 2 3 1
1 2 3 1
1 2 3 1
"""
parse(PHononInput,content)

end