# LTR Consensus Sequence Builder for RepeatMasker

## Overview
This toolkit converts Look4LTRs RTR output files into RepeatMasker-compatible consensus sequence libraries for genome annotation.

## RTR File Format
The RTR (RetroTransposon Report) file from Look4LTRs contains the following columns:
- **chrom**: Chromosome/contig name
- **ID**: Element ID
- **LeftStart/LeftEnd**: 5' LTR coordinates
- **RightStart/RightEnd**: 3' LTR coordinates
- **GraphGroup**: Family clustering group
- **LTRIdentity**: Identity between 5' and 3' LTRs
- **CaseType**: Element type (Single, SoloSingle, etc.)
- Additional fields: PPT, TSD coordinates, nesting information

## Available Scripts

### 1. Quick Bash Pipeline: `rtr2repeatmasker.sh`
**Fastest and simplest option for basic consensus building**

```bash
./rtr2repeatmasker.sh genome.rtr genome.fa output_dir/ 8
```

Features:
- Extracts LTR sequences from genome
- Groups by GraphGroup
- Aligns with MAFFT/MUSCLE
- Builds simple majority-rule consensus
- Creates RepeatMasker library file

Output:
- `LTR_library.lib`: RepeatMasker library
- `library_stats.txt`: Family statistics
- `run_repeatmasker.sh`: Ready-to-use RepeatMasker script

### 2. Python Basic: `build_ltr_consensus.py`
**More control over parameters and filtering**

```bash
python3 build_ltr_consensus.py \
    -r genome.rtr \
    -g genome.fa \
    -o output_dir/ \
    --min-identity 0.8 \
    --min-elements 3 \
    --threads 8
```

Features:
- Filters by LTR identity threshold
- Extracts complete elements with flanking regions
- Multiple sequence alignment
- Consensus with ambiguity codes
- Detailed statistics output

### 3. Python Advanced: `build_ltr_consensus_advanced.py`
**More comprehensive with CD-HIT clustering**

```bash
python3 build_ltr_consensus_advanced.py \
    -r genome.rtr \
    -g genome.fa \
    -o output_dir/ \
    --min-identity 0.8 \
    --min-elements 3 \
    --cluster-identity 0.9
```

Features:
- CD-HIT clustering for better family definition
- Handles reverse complement elements
- Extracts TSD and PPT regions
- Advanced consensus building
- Multiple alignment strategies

### 4. Python Enhanced: `build_ltr_consensus_enhanced.py` ⭐
**Most sophisticated with full RTR information usage**

```bash
python3 build_ltr_consensus_enhanced.py \
    -r genome.rtr \
    -g genome.fa \
    -o output_dir/ \
    --min-identity 0.7 \
    --min-elements 2
```

**New Enhanced Features:**
- **LTRIdentity-weighted consensus**: Higher identity sequences have more influence
- **NestIn/NestOut relationship parsing**: Handles nested element structures properly
- **Quality-tier classification**: Excellent (≥0.95), High (≥0.9), Medium (≥0.8), Low (≥0.7)
- **Subfamily structure**: Creates subfamilies based on identity clustering
- **CaseType-aware processing**: Different handling for Single, SoloSingle, Nested elements
- **Canonical library**: High-quality subset for stringent annotation
- **Enhanced statistics**: Detailed subfamily composition and quality metrics

## Workflow Pipeline

```
RTR File + Genome FASTA
         ↓
    Extract Sequences
         ↓
    Group by Family
         ↓
    Align Sequences
         ↓
   Build Consensus
         ↓
  RepeatMasker Library
```

## Using the Library with RepeatMasker

### Basic Usage
```bash
RepeatMasker -lib LTR_library.lib -pa 8 genome.fasta
```

### Advanced Usage
```bash
RepeatMasker \
    -lib LTR_library.lib \
    -pa 8 \
    -xsmall \        # Soft-mask (lowercase) repeats
    -html \          # Create HTML output
    -gff \           # Create GFF3 output
    -dir output/ \   # Output directory
    genome.fasta
```

### Combine with RepBase
```bash
# Concatenate with RepBase library
cat LTR_library.lib /path/to/RepBase.lib > combined.lib
RepeatMasker -lib combined.lib genome.fasta
```

## Requirements

### Essential
- Python 3.6+ with BioPython
- samtools (for sequence extraction)
- MAFFT or MUSCLE (for alignment)

### Optional
- CD-HIT (for clustering)
- pandas, numpy (for advanced script)
- RepeatMasker (for annotation)

### Installation
```bash
# Install Python dependencies
pip install biopython pandas numpy

# Install alignment tools (choose one)
conda install -c bioconda mafft
# OR
conda install -c bioconda muscle

# Optional: Install CD-HIT
conda install -c bioconda cd-hit
```

## Enhanced RTR Information Usage

The enhanced script (`build_ltr_consensus_enhanced.py`) fully utilizes all RTR file columns:

### **LTRIdentity Usage**
- **Quality weighting**: Elements with higher LTRIdentity contribute more to consensus
- **Subfamily creation**: Elements grouped by identity tiers (Excellent ≥0.95, High ≥0.9, etc.)
- **Library separation**: Canonical library contains only high-quality elements

### **NestIn/NestOut Processing**
- **Nested relationship parsing**: Understands which elements are nested within others
- **Host element identification**: Elements with NestOut are treated as hosts
- **Clean internal regions**: Only uses non-nested elements for internal consensus
- **Nested element tracking**: Statistics include nesting complexity

### **CaseType-Aware Processing**
- **Single**: Complete elements used for full consensus sequences
- **SoloSingle**: Only LTRs used, no internal region
- **RecentlyNested**: Special handling for recently nested elements
- **Different extraction**: Each type processed according to its structure

### **Additional Features**
- **TSD extraction**: Target Site Duplications when available
- **PPT extraction**: Polypurine Tracts when available  
- **Reverse complement handling**: Proper orientation for RC elements
- **GraphGroup optimization**: Enhanced family definition beyond simple grouping

## Output Files

### Enhanced Library Files
```
# Main library with all subfamilies
>FAM22_high_0_LTR#LTR/Gypsy elements=15 avg_id=0.921
TGCAGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT

# Canonical library (high quality only)  
>FAM22_excellent_0_LTR#LTR/Gypsy elements=8 avg_id=0.965
TGCAGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
```

### Enhanced Statistics
```
Subfamily        GraphGroup  NumElements  ExcellentElements  HighElements  QualityTiers
FAM22_high_0     22          15           5                  10            {'excellent':5,'high':10}
FAM29_medium_0   29          8            0                  3             {'high':3,'medium':5}
```

### Nesting Analysis
```
Total elements: 1,234
Nested elements: 156 (12.6%)
Host elements: 89 (7.2%)
Complex nesting (depth >2): 23 (1.9%)
```

## Tips and Best Practices

1. **Filter by Identity**: Use `--min-identity 0.8` to include only high-quality elements
2. **Minimum Elements**: Set `--min-elements 3` to ensure reliable consensus
3. **Check Alignment Tools**: Ensure MAFFT or MUSCLE is installed before running
4. **Memory Usage**: For large genomes, process chromosomes separately
5. **Validation**: Check consensus sequences with multiple alignment viewer

## Troubleshooting

### No consensus sequences generated
- Check if RTR file has sufficient high-identity elements
- Lower `--min-identity` threshold
- Ensure genome FASTA headers match RTR chromosome names

### Alignment fails
- Install MAFFT: `conda install -c bioconda mafft`
- Check sequence lengths (very short sequences may fail)

### CD-HIT not working
- Install: `conda install -c bioconda cd-hit`
- Use basic script without clustering

## Example Complete Workflow

```bash
# 1. Run Look4LTRs
./look4ltrs -f genome.fa -o look4ltrs_output/ -pa 8

# 2. Build consensus library
./rtr2repeatmasker.sh \
    look4ltrs_output/genome.rtr \
    genome.fa \
    repeatmasker_lib/ \
    8

# 3. Run RepeatMasker
cd repeatmasker_lib/
./run_repeatmasker.sh ../genome.fa 8

# 4. View results
less RM_output/genome.fa.out
```

## Citation
If you use these scripts, please cite:
- Look4LTRs: [Original publication]
- RepeatMasker: Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0.

## Support
For issues or questions:
- Look4LTRs: https://github.com/[repository]
- Scripts: Created by Assistant for LTR consensus building

---
*Last updated: 2025-01-13*