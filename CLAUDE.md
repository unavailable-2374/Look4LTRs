# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Look4LTRs is a C++ tool for detecting Long Terminal Repeat Retrotransposons (LTR RTs) in genomic sequences. It can find recently nested LTR RTs de novo and is designed for both model and non-model organisms.

## Build System

The project uses CMake with the following requirements:
- GNU g++ 7.5.0 or later (11.1.0 recommended)
- cmake 3.10.3 or later
- C++17 standard

### Core Build Commands

```bash
# Initial setup
mkdir bin
cd bin

# Configure build
cmake ..
# Or with specific compiler:
cmake .. -DCMAKE_CXX_COMPILER=g++-11

# Build main executable
make look4ltrs

# Build training executables
make generateTrainingData
make generateGraphData

# Other useful executables
make red          # RED repeat detector
make identity     # Identity calculator
make meshclust    # Clustering tool
```

## Code Architecture

### Core Modules

1. **LTR Detection Pipeline** (`src/ltr/`):
   - `look4ltrs.cpp` - Main entry point
   - `ModulePipeline.h/cpp` - Orchestrates the detection pipeline
   - `Detector.h/cpp` - Core detection logic using SGD classifier
   - `Matcher.h/cpp` - LTR matching algorithms
   - `Element.h/cpp` - Represents detected LTR elements

2. **Machine Learning Components** (`src/`):
   - `SGD.h/cpp` - Stochastic Gradient Descent classifier
   - `GLM*.h/cpp` - Generalized Linear Model components
   - `FeatureExpander.h/cpp` - Feature engineering
   - `Normalizer.h/cpp` - Data normalization

3. **Output Formats** (`src/ltr/Output*.h/cpp`):
   - `OutputBed.h/cpp` - BED format output
   - `OutputRtr.h/cpp` - RTR (RetroTransposon Relationship) format
   - `OutputCpx.h/cpp` - Complex regions format

4. **Utility Libraries**:
   - `src/utility/` - Location handling and utilities
   - `src/exception/` - Custom exception classes
   - `src/meshclust/` - Clustering algorithms
   - `src/nonltr/` - Non-LTR repeat detection (RED)

### Key Classes and Parameters

- `LtrParameters` - Contains all detection parameters (K=13, distance thresholds)
- `DirectedGraph` - Graph-based LTR family grouping
- `DeepNesting` - Handles nested LTR detection
- `FastaReader` - FASTA file parsing

## Usage Patterns

### Basic LTR Detection
```bash
./look4ltrs --fasta /path/to/genome.fa --out /output/dir/ --parallel 8
```

### Training Mode
```bash
./look4ltrs --fasta /predict/genome/ --train /training/genome/ --out /output/ --parallel 8
```

### Custom Model
```bash
./look4ltrs --fasta /genome/ --out /output/ --config /custom/config.txt
```

## Training Pipeline

Python training components are in `Training/`:
- `trainModel.py` - Main training script
- `requirements.txt` - Python dependencies (biopython, scikit-learn, etc.)
- Uses SGD classifier from scikit-learn
- Requires BED format annotations for training data

## Output Analysis Scripts

Scripts in `scripts/` for result analysis:
- `findRecentNest.py` - Find recently nested LTR RTs
- `findSameGraphNest.py` - Find same-family nested elements  
- `findRT.py` - Look up specific LTR RT by ID

## File Requirements

- Input: FASTA files with `.fa` extension only
- Multi-FASTA supported but memory-intensive
- Output creates three directories: Bed/, Rtr/, Cpx/

## Testing

The project includes test executables but no automated test framework:
- Test files in `src/test/`
- Build tests with `make chromtest`, `make driver`, etc.
- Manual testing approach - check for specific test commands in build logs

## Development Notes

- Heavy use of forward declarations (`.fwd.h` files)
- OpenMP for parallelization
- Exception-based error handling
- Static libraries for modular compilation
- No linting/formatting tools configured - follow existing code style