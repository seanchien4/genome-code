# braker_to_chrom.py
gtf_file = "braker.gtf"
chrom_out = "genome.chrom"

# Parse GTF, get one entry per transcript with min(start) and max(end) across CDS
genes = {}  # transcript_id -> [scaffold, strand, start, stop]

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        scaffold, source, feature, start, end, score, strand, frame, attr = fields
        if feature != "CDS":
            continue
        # extract transcript_id
        tid = None
        for a in attr.split(";"):
            a = a.strip()
            if a.startswith("transcript_id"):
                tid = a.split('"')[1] if '"' in a else a.split()[1]
                break
        if tid is None:
            continue
        start, end = int(start), int(end)
        if tid not in genes:
            genes[tid] = [scaffold, strand, start, end]
        else:
            genes[tid][2] = min(genes[tid][2], start)
            genes[tid][3] = max(genes[tid][3], end)

with open(chrom_out, "w") as out:
    for tid, (scaf, strand, start, stop) in genes.items():
        out.write(f"{tid}\t{scaf}\t{strand}\t{start}\t{stop}\n")

print(f"Wrote {len(genes)} entries to {chrom_out}")