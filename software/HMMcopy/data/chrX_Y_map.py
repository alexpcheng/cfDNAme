import statistics
with open('500.hg19.full.map.wig') as f:
    chrX_maps = []
    chrY_maps = []
    chrY=False
    chrX=False
    for line in f:
        value = line.strip()
        if value == 'fixedStep chrom=chrX start=1 step=500 span=500':
            chrX=True
            chrY=False
            continue
        if value == 'fixedStep chrom=chrY start=1 step=500 span=500':
            chrX=False
            chrY=True
            continue
        if chrX:
            chrX_maps.append(float(value))
        if chrY:
            chrY_maps.append(float(value))
mean_map_X = statistics.mean(chrX_maps)
mean_map_Y = statistics.mean(chrY_maps)
print(mean_map_X)
print(mean_map_Y)
