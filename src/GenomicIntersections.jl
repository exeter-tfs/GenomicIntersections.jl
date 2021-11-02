module GenomicIntersections


using DataFrames, GenomicFeatures, Statistics

export countoverlap, intersecttable, overlappinglocations, peakintersect, markintersection!, annotatecol!

include("intersections.jl")



end
