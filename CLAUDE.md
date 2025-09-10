# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Look4LTRs is a Long Terminal Repeat (LTR) Retrotransposon detection tool capable of finding recently nested LTR RTs de novo. It's written in C++17 with Python scripts for analysis and training.

## Build System

The project uses CMake (minimum version 3.10) with specific compiler requirements:
- GNU g++ 11.1.0 or later (7.5.0 minimum)
- C++17 standard required
- OpenMP enabled for parallelization

### Common Build Commands

```bash
# Initial setup
mkdir bin
cd bin

# Standard build
cmake ..
make look4ltrs

# With specific compiler
cmake .. -DCMAKE_CXX_COMPILER=g++-11

# Build training executables
make generateTrainingData
make generateGraphData

# Other useful executables
make red              # RED repeat detector
make meshclust        # Clustering utility
make identity         # Identity calculator
```

## Architecture Overview

The codebase is organized into several key modules:

### Core Libraries (Static)
- **main**: Core functionality including FASTA reading, data generation, GLM classifiers, feature processing
- **ltr**: LTR-specific detection pipeline including matching, filtering, output generation
- **clustering**: MeShClust clustering algorithms and evaluation
- **nonltr**: Non-LTR repeat detection (RED integration)
- **utility**: Common utilities and location handling
- **exception**: Custom exception classes

### Key Components

1. **ModulePipeline** (`src/ltr/ModulePipeline.{h,cpp}`): Main orchestration class that coordinates the LTR detection pipeline
2. **Detector** (`src/ltr/Detector.{h,cpp}`): Machine learning-based element detection using SGD classifier
3. **Matcher** (`src/ltr/Matcher.{h,cpp}`): LTR matching and alignment logic
4. **RED Integration** (`src/red/Red.{h,cpp}`): Non-LTR repeat detection system
5. **Output Classes**: BED, RTR, and CPX format writers in `src/ltr/Output*.{h,cpp}`

### Processing Flow
1. FASTA input processing
2. K-mer histogram analysis
3. RED preprocessing for repeat detection
4. LTR stretch identification and scoring
5. Element matching and filtering
6. Nested element detection (DeepNesting)
7. Output generation in multiple formats

## Python Components

### Analysis Scripts (`scripts/`)
- `findRecentNest.py`: Find recently nested LTR RTs
- `findSameGraphNest.py`: Find same-family nested elements  
- `findRT.py`: Lookup specific RT elements by ID

### Training Pipeline (`Training/`)
- `trainModel.py`: Main training script for custom models
- `generateSemiSynthetic.py`: Generate synthetic training data
- `detectorTrainer.py`: SGD classifier training
- Requirements: biopython, numpy, scikit-learn, scipy

## Usage Patterns

### Standard Detection
```bash
./look4ltrs --fasta /path/to/genome.fa --out /output/dir --parallel 8
```

### Multi-genome Analysis
```bash
./look4ltrs --fasta /genome1/dir/ /genome2/dir/ --out /output/ --parallel 8
```

### Training Mode
```bash
./look4ltrs --fasta /test/genome/ --train /training/genome/ --out /output/ --parallel 8
```

## File Requirements

- Input FASTA files MUST have `.fa` extension
- Training requires BED files with 7 columns: chrom, start, end, left_start, left_end, right_start, right_end
- Multi-FASTA supported but memory-intensive

## Output Formats

- **BED**: Standard 3-column format (chrom, start, end)
- **RTR**: 16-column Look4LTRs format with detailed element information
- **CPX**: Variable-column format for complex/uncertain regions

## Memory and Performance

- Designed for whole-genome analysis
- OpenMP parallelization available (`--parallel` parameter)
- Memory usage scales with genome size and thread count
- Consider reducing threads if memory constraints occur

## Known Issues and Fixes

### Deep Nesting Infinite Loop
- **Problem**: Software may get stuck in "Found a nest at level XXX" loops due to infinite recursion in `DeepNesting::findDeep()`
- **Location**: `src/ltr/DeepNesting.cpp:127`
- **Fix Applied**: 
  - Added maximum recursion depth limit as class static constant (MAX_NESTING_DEPTH = 10)
  - Added visited regions tracking to prevent processing same region multiple times
  - Fixed OpenMP critical block compatibility issues
- **Symptoms**: Console output continuously shows increasing nest levels without termination

### Multithreading Not Working
- **Problem**: Program runs single-threaded even when `--parallel` is specified with multiple threads
- **Original Issue**: Parallelization was at FASTA file level, not chromosome level
- **Fix Applied**: Modified `src/ltr/look4ltrs.cpp` to parallelize at chromosome level
- **Details**: 
  - Collects all chromosomes from all FASTA files first
  - Then processes chromosomes in parallel using OpenMP
  - Shows message "Processing X chromosomes with Y threads..."