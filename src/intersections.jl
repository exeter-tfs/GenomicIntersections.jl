

intervals(table; off = 0) = @with table IntervalCollection(Interval.(table.chrom, table.start .+ off, table.stop, '.', 1:size(table, 1)))

@inline unionlocation(a, b) = min(a.start, b.start):max(a.stop, b.stop)
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