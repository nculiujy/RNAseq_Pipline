import re, sys

gtf, bed = sys.argv[1], sys.argv[2]
transcripts, exons = {}, {}

with open(gtf) as f:
    for line in f:
        if line.startswith("#"): continue
        fields = line.strip().split("\t")
        if len(fields) < 9: continue
        chrom, _, ftype, start, end, _, strand, _, attrs = fields
        start, end = int(start)-1, int(end)
        tid = re.search(r'transcript_id "([^"]+)"', attrs)
        if not tid: continue
        tid = tid.group(1)
        if ftype == "transcript":
            transcripts[tid] = [chrom, start, end, tid, 0, strand]
        elif ftype == "exon":
            exons.setdefault(tid, []).append((start, end))

with open(bed, "w") as out:
    for tid, (chrom, tstart, tend, name, score, strand) in transcripts.items():
        ex = sorted(exons.get(tid, [(tstart, tend)]))
        sizes = ",".join(str(e-s) for s, e in ex) + ","
        starts = ",".join(str(s-tstart) for s, e in ex) + ","
        out.write(f"{chrom}\t{tstart}\t{tend}\t{name}\t{score}\t{strand}\t{tstart}\t{tend}\t0\t{len(ex)}\t{sizes}\t{starts}\n")

print(f"Done: {len(transcripts)} transcripts")
