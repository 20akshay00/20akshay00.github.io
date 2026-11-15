# Shared look for the plots in the tidbits, included from each post's (hidden) setup cell.
# The page text is Literata, which GR cannot load; Palatino is its closest built-in serif.
# The sizes are a little larger than the defaults because the figures are scaled down to
# fit the notebook page, which would otherwise make the tick labels tiny.
default(
    fontfamily = "palatino",
    titlefontsize = 13,
    guidefontsize = 12,
    tickfontsize = 10,
    legendfontsize = 10,
)
