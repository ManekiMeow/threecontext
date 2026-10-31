import glob, re, statistics, os
SC = "/tmp/claude-0/-home-user-threecontext/5e407493-7848-5a52-a504-a0152de6dc76/scratchpad/night"
for f in sorted(glob.glob(SC + "/*_cayley4_candidate_*.log")) + sorted(glob.glob(SC + "/*control*.log")):
    vals = []
    for line in open(f):
        m = re.search(r"\|\|G\^2 - 4G\|\| = ([0-9.e+-]+)", line) or re.search(r"objective ([0-9.e+-]+)", line)
        if m and "restart" in line:
            vals.append(float(m.group(1)))
    if vals:
        print(f"{os.path.basename(f):60s} starts {len(vals):4d}  best {min(vals):.4g}  median {statistics.median(vals):.4g}  worst {max(vals):.4g}")
