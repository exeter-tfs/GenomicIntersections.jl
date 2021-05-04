

intervals(table; off = 0) = IntervalCollection(Interval.(table.chrom, table.start .+ off, table.stop, '.', 1:size(table, 1)))
span(interval) = rightposition(interval) - leftposition(interval) + 1
intersectinginterval(ia, ib) = Interval(seqname(ia), max(leftposition(ia), leftposition(ib)), min(rightposition(ia), rightposition(ib)))
@inline unionlocation(a, b) = min(a.start, b.start):max(a.stop, b.stop)


function countoverlap(tableA, tableB)
    ivA = intervals(tableA)::IntervalCollection{Int64}
    ivB = intervals(tableB)::IntervalCollection{Int64}

    spanA = sum(span, ivA)
    spanB = sum(span, ivB)
    spanI = 0
    n = 0
    totalA = Set{Int64}()
    totalB = Set{Int64}()
    for (ia, ib) in eachoverlap(ivA, ivB)
        ii = intersectinginterval(ia, ib)
        spanI += span(ii)
        push!(totalA, metadata(ia))
        push!(totalB, metadata(ib))
        n += 1
    end

    jaccard = spanI/(spanA + spanB - spanI)

    n, length(totalA), length(totalB), jaccard
end

function intersecttable(tableA, tableB)
    ivA = intervals(tableA)::IntervalCollection{Int64}
    ivB = intervals(tableB)::IntervalCollection{Int64}
    itable = DataFrame(chrom=String[], start=Int[], stop=Int[], name=String[], score=Int[], strand=String[], FC=Float64[], nlogp=Float64[], nlogq=Float64[], summit=Int[])
    for (ia, ib) in eachoverlap(ivA, ivB)
        rA = tableA[metadata(ia), :]
        rB = tableB[metadata(ib), :]

        push!(itable, (rA.chrom, max(rA.start, rB.start), min(rA.stop, rB.stop), "", div(rA.score + rB.score, 2), ".", (rA.FC + rB.FC)/2, max(rA.FC, rB.FC), max(rA.FC, rB.FC), div(rA.summit + rB.summit, 2)))
    end
    itable
end


function overlappinglocations(chroms, locations)

    regions = zeros(Int, length(chroms))
    regionindex = 0
    location = 1:0
    chrom    = ""

    for i = 1:length(locations)
        if (chroms[i] == chrom) && !isempty(intersect(location, locations[i]))
            regions[i] = regionindex
            location = unionlocation(location, locations[i])
        else
            regionindex += 1
            regions[i] = regionindex
            location = locations[i]
            chrom  = chroms[i]
        end
    end
    regions
end


peakintersect(labels::Vector{T}, peakdict::Dict{V, DataFrame}) where {T, V} = peakintersect(labels, [peakdict[l] for l in labels])
function peakintersect(labels::Vector{T}, peaks::Vector{DataFrame}) where {T}
    allpeaks = reduce(vcat, peaks)

    sort!(allpeaks, [:chrom, :start, :stop, :Origin]);
    chroms    = allpeaks.chrom::Vector{String}
    locations = [(a+1):b for (a, b) in zip(allpeaks.start, allpeaks.stop)]::Vector{UnitRange{Int}}
    allpeaks.Group = overlappinglocations(chroms, locations);


    combpeaks = combine(groupby(allpeaks, :Group),
                    :chrom => first => :chrom,
                    :start => minimum => :start,
                    :stop => minimum => :stop,
                    :nrow => :TotalPeaks,
                    :score => mean => :score,
                    :Origin => Ref => :Origin)

    membership =  DataFrame([1(l .∈ combpeaks.Origin) for l in labels], Symbol.(labels))

    [comb_peaks[!, Not([:Group, :Origin])] membership]
end