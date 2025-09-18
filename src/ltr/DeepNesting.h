/*
 * LtrDetector v2.0 annotates LTR retro-transposons in a genome.
 * 
 * DeepNesting
 * 
 *  Created on: X X, 20XX
 *      Author: Anthony B. Garza.
 * Reviewer:
 *   Purpose: 
 *        
 * 
 * Academic use: Affero General Public License version 1.
 *
 * Any restrictions to use for-profit or non-academics: Alternative commercial license is needed.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * Please contact Dr. Hani Z. Girgis (hzgirgis@buffalo.edu) if you need more information.
 * 
 * Copyright (C) 2022 by the authors.
 */

#pragma once

#include "ModulePipeline.fwd.h"
#include "DeepNesting.fwd.h"

#include "ModulePipeline.h"
#include "RT.h"
#include "../red/Red.h"
#include "../IdentityCalculator.h"

#include <vector>
#include <string>
#include <set>
#include <utility>


class DeepNesting
{
private:
    // Variables
    std::vector<RT *> *rtVecPtr;
    const std::string *chrom;
    ModulePipeline &mp;
    std::vector<ModulePipeline*> mpVec;
    std::vector<RT*> recentNestVec;
    Red &red;
    IdentityCalculator<int32_t> &icStandard;
    IdentityCalculator<int32_t> &icRecent;
    
    // Maximum nesting depth limit to prevent infinite recursion
    static constexpr int MAX_NESTING_DEPTH = 3;  // Reduced to 3 for faster processing
    
    // Track visited regions to prevent infinite loops with region hash
    std::set<std::pair<int, int>> visitedRegions;
    
    // Cache for processed regions to avoid redundant work
    std::map<std::pair<int, int>, std::vector<RT*>> regionCache;
    
    // Early termination threshold for complex regions
    static constexpr int MAX_SAME_REGION_VISITS = 2;  // Reduced to 2 to prevent excessive loops

    // Methods

public:
    
    // Constructor

    DeepNesting(std::vector<RT*>* _rtVecPtr, const std::string *_chrom, ModulePipeline &_mp, Red &_red, IdentityCalculator<int32_t> &_icStandard, IdentityCalculator<int32_t> &_icRecent);
    ~DeepNesting();

    // Methods
    void findRegion();
    std::vector<ModulePipeline*> findDeep(std::string &graphSeq, int removeStart, int removeLength, int level);
    
    // Helper method to check if region should be skipped
    bool shouldSkipRegion(int removeStart, int removeLength, int level);
    
    // Track how many times each region has been visited
    std::map<std::pair<int, int>, int> regionVisitCount;
    
    // Global counter for total regions processed to prevent runaway recursion
    int totalRegionsProcessed = 0;
    static constexpr int MAX_TOTAL_REGIONS = 50;  // Maximum total regions to process
};