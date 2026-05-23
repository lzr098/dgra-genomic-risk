import re

base = '/root/.openclaw/skills/dgra-genomic-risk'

with open(f'{base}/scripts/dgra_core.py', 'r') as f:
    content = f.read()

# Resolve merge conflicts by keeping HEAD (remote) features + our fixes
# Conflict 1: around line 4452 - gnomAD batch query section
# HEAD has MyVariant.info integration; we need to keep it

# Strategy: replace all conflict markers with HEAD content,
# then manually re-apply our NMD fix

# Remove all conflict markers - keep HEAD content between <<<<<<< HEAD and =======,
# discard >>>>>>> branch content

def resolve_conflicts(text):
    """Keep HEAD (remote) version of all conflicts."""
    pattern = r'<<<<<<< HEAD\n(.*?)=======\n.*?>>>>>>> .*?\n'
    
    def repl(match):
        head_content = match.group(1)
        return head_content
    
    return re.sub(pattern, repl, text, flags=re.DOTALL)

resolved = resolve_conflicts(content)

# Now re-apply our NMD fix: find the NMD section and add gene_constraint write
# Look for: v.nmd_prediction = predict_nmd(...)
nmd_pattern = r'(v\.nmd_prediction = predict_nmd\([^)]+\))\n\s+nmd_count \+= 1'
nmd_replacement = r'''\1
            # v0.9.3: Write NMD prediction into gene_constraint for JSON report access
            if v.gene_constraint is None:
                v.gene_constraint = {}
            v.gene_constraint["nmd_prediction"] = v.nmd_prediction
            nmd_count += 1'''

resolved = re.sub(nmd_pattern, nmd_replacement, resolved)

with open(f'{base}/scripts/dgra_core.py', 'w') as f:
    f.write(resolved)

# Verify no remaining conflict markers
remaining = re.findall(r'<<<<<<<|=======|>>>>>>>', resolved)
if remaining:
    print(f"WARNING: Remaining conflict markers: {len(remaining)}")
else:
    print("All conflicts resolved")
