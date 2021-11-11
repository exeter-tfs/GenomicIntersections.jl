module GenomicIntersections


using DataFrames, GenomicFeatures, Statistics

export countoverlap, intersecttable, overlappinglocations, peakintersect, markintersection!, annotatecol!, countintersection!

include("intersections.jl")



end
