module PHonon

using PyFortran90Namelists: Parser
using ReadableRegex: RegexString, char_in, exactly, one_or_more, zero_or_more, capture

export parse_frequency, parse_dos

const AD = "[0-9]"
const NEWLINE = RegexString(raw"\R")
const INTEGER = RegexString("([-+]?$AD+)")
const FIXED_POINT_REAL = RegexString("([-+]?$AD*\\.$AD+|$AD+\\.?$AD*)")
const FREQ_PER_ROW =
    exactly(6, one_or_more(char_in(" \t")) * FIXED_POINT_REAL) * zero_or_more(NEWLINE)
const QCOORD = capture(
    exactly(3, one_or_more(char_in(" \t")) * FIXED_POINT_REAL) * one_or_more(NEWLINE);
    as = "coord",
)
_QBOLCK(div, rem) =
    QCOORD * capture(
        exactly(div, FREQ_PER_ROW) *
        exactly(rem, one_or_more(char_in(" \t")) * FIXED_POINT_REAL) *
        zero_or_more(NEWLINE);
        as = "band",
    )

function parse_frequency(str::AbstractString)
    d::Dict{String,Any} = Parser().reads(str)["plot"]
    nks, nbnd = d["nks"], d["nbnd"]
    nrows, remcols = divrem(nbnd, 6)  # QE splits branches into 6 columns per line
    regex = _QBOLCK(nrows, remcols)
    vec = map(eachmatch(regex, str)) do m
        coord, band = m[:coord], m[:band]
        fcoord = Tuple(parse(Float64, x) for x in split(coord, " "; keepempty = false))
        fband = Tuple(parse(Float64, y) for y in split(band, " "; keepempty = false))
        fcoord => fband
    end
    @assert length(vec) == nks
    return vec
end

function parse_dos(str::AbstractString)
    map(split(str, r"\R"; keepempty = false)) do line
        x, dos = [parse(Float64, x) for x in split(line, " "; keepempty = false)]
    end
end

end
