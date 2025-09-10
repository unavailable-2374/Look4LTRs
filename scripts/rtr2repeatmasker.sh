#!/bin/bash

# Quick RTR to RepeatMasker library converter
# Usage: ./rtr2repeatmasker.sh <rtr_file> <genome_fasta> <output_dir>

set -e

RTR_FILE=$1
GENOME_FASTA=$2
OUTPUT_DIR=$3
THREADS=${4:-8}

if [ $# -lt 3 ]; then
    echo "Usage: $0 <rtr_file> <genome_fasta> <output_dir> [threads]"
    echo "Example: $0 genome.rtr genome.fa output/ 8"
    exit 1
fi

echo "=========================================="
echo "RTR to RepeatMasker Library Converter"
echo "=========================================="
echo "RTR file: $RTR_FILE"
echo "Genome: $GENOME_FASTA"
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "=========================================="

# Create output directory
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# Step 1: Extract LTR sequences based on RTR coordinates
echo "[$(date)] Step 1: Extracting LTR sequences..."

awk 'NR>1 {
    # Skip header
    if ($16 >= 0.7) {  # Filter by LTR identity
        print $1":"$3"-"$4" "$1"_"$2"_L "$15" "$16
        print $1":"$5"-"$6" "$1"_"$2"_R "$15" "$16
    }
}' $RTR_FILE | while read coords id group identity; do
    echo ">$id group=$group identity=$identity"
    samtools faidx $GENOME_FASTA $coords | grep -v ">"
done > ltr_sequences.fasta

echo "[$(date)] Extracted $(grep -c ">" ltr_sequences.fasta) LTR sequences"

# Step 2: Cluster sequences by GraphGroup
echo "[$(date)] Step 2: Clustering by GraphGroup..."

mkdir -p groups
awk '/^>/ {
    group = $2
    gsub("group=", "", group)
    file = "groups/group_" group ".fasta"
    print > file
    getline
    print > file
}' ltr_sequences.fasta

# Step 3: Build consensus for each group
echo "[$(date)] Step 3: Building consensus sequences..."

> consensus_library.fasta

for group_file in groups/group_*.fasta; do
    if [ ! -s "$group_file" ]; then
        continue
    fi
    
    group=$(basename $group_file .fasta | sed 's/group_//')
    num_seqs=$(grep -c ">" $group_file)
    
    if [ $num_seqs -lt 4 ]; then
        echo "  Skipping group $group (only $num_seqs sequences)"
        continue
    fi
    
    echo "  Processing group $group ($num_seqs sequences)..."
    
    # Align sequences
    if command -v mafft &> /dev/null; then
        mafft --thread $THREADS --quiet $group_file > groups/group_${group}_aligned.fasta 2>/dev/null
    elif command -v muscle &> /dev/null; then
        muscle -in $group_file -out groups/group_${group}_aligned.fasta -quiet 2>/dev/null
    else
        echo "Warning: No alignment tool found (mafft or muscle required)"
        continue
    fi
    
    # Build consensus using simple majority rule
    if [ -s "groups/group_${group}_aligned.fasta" ]; then
        python3 - <<EOF >> consensus_library.fasta
import sys
from collections import Counter

def build_consensus(aligned_file):
    sequences = []
    with open(aligned_file) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    
    if not sequences:
        return ""
    
    consensus = []
    for i in range(len(sequences[0])):
        column = [seq[i] for seq in sequences if i < len(seq)]
        bases = [b for b in column if b != '-' and b != 'N']
        if bases:
            most_common = Counter(bases).most_common(1)[0][0]
            consensus.append(most_common)
    
    return ''.join(consensus)

consensus = build_consensus("groups/group_${group}_aligned.fasta")
if consensus:
    print(">LTR_FAM${group} # LTR/Unknown")
    for i in range(0, len(consensus), 60):
        print(consensus[i:i+60])
EOF
    fi
done

echo "[$(date)] Created $(grep -c ">" consensus_library.fasta) consensus sequences"

# Step 4: Optional CD-HIT clustering for redundancy removal
if command -v cd-hit-est &> /dev/null; then
    echo "[$(date)] Step 4: Removing redundancy with CD-HIT..."
    cd-hit-est -i consensus_library.fasta -o consensus_library_nr.fasta \
        -c 0.9 -aS 0.8 -T $THREADS -M 2000 -d 0 &> cdhit.log
    mv consensus_library_nr.fasta consensus_library.fasta
    echo "[$(date)] Redundancy removal complete"
else
    echo "[$(date)] Step 4: Skipping CD-HIT (not installed)"
fi

# Step 5: Create RepeatMasker-compatible library
echo "[$(date)] Step 5: Creating RepeatMasker library..."

# Add RepeatMasker classification tags
sed 's/# LTR\/Unknown/#LTR/' consensus_library.fasta > LTR_library.lib

# Create simple statistics
echo "Library Statistics" > library_stats.txt
echo "==================" >> library_stats.txt
echo "Total families: $(grep -c ">" LTR_library.lib)" >> library_stats.txt
echo "" >> library_stats.txt
echo "Family ID       Length" >> library_stats.txt
echo "--------        ------" >> library_stats.txt
awk '/^>/ {
    if (seq) print name, length(seq)
    name = substr($1, 2)
    seq = ""
    next
}
{seq = seq $0}
END {if (seq) print name, length(seq)}' LTR_library.lib >> library_stats.txt

# Create RepeatMasker usage script
cat > run_repeatmasker.sh <<'SCRIPT'
#!/bin/bash

# RepeatMasker usage with custom LTR library

GENOME=$1
THREADS=${2:-8}

if [ -z "$GENOME" ]; then
    echo "Usage: $0 <genome.fasta> [threads]"
    exit 1
fi

echo "Running RepeatMasker with custom LTR library..."
RepeatMasker \
    -lib LTR_library.lib \
    -pa $THREADS \
    -xsmall \
    -html \
    -gff \
    -dir RM_output \
    $GENOME

echo "RepeatMasker complete. Results in RM_output/"
SCRIPT

chmod +x run_repeatmasker.sh

# Final summary
echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo "Output files:"
echo "  - LTR_library.lib: RepeatMasker library with $(grep -c ">" LTR_library.lib) families"
echo "  - library_stats.txt: Family statistics"
echo "  - run_repeatmasker.sh: Script to run RepeatMasker"
echo ""
echo "To annotate a genome:"
echo "  cd $OUTPUT_DIR"
echo "  ./run_repeatmasker.sh your_genome.fasta"
echo "=========================================="