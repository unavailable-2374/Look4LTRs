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

#include "DeepNesting.h"
#include <map>

DeepNesting::DeepNesting(std::vector<RT*>* _rtVecPtr, const std::string *_chrom, ModulePipeline &_mp, Red &_red, IdentityCalculator<int32_t> &_icStandard, IdentityCalculator<int32_t> &_icRecent) : rtVecPtr(_rtVecPtr), chrom(_chrom), mp(_mp), red(_red), icStandard(_icStandard), icRecent(_icRecent)
{

}

DeepNesting::~DeepNesting()
{
    for (auto &mp : mpVec) {
        delete mp;
        mp = nullptr;
    }
}

void DeepNesting::findRegion() {
    // Reset counters for each new region search
    visitedRegions.clear();
    regionCache.clear();
    totalRegionsProcessed = 0;
    
    for (auto rt : *rtVecPtr) {
        // Searching for the recently nested inside element
        if (rt->getCaseType() == "RecentlyNestedInner") {

            // Get sequence of the graph that these recently nested elements are in
            std::pair<int, int> graphRegion = mp.getFamilyRegion(rt); // start and end of graph
            std::string graphSeq = chrom->substr(graphRegion.first, graphRegion.second - graphRegion.first);

            // Where the nested element starts in the graph
            int nestStart = rt->getStart() - graphRegion.first;
            // How long is the nested element?
            int length = rt->getSize();

            auto nestModuleVec = findDeep(graphSeq, nestStart, length, 2);
            mpVec.insert(mpVec.end(), nestModuleVec.begin(), nestModuleVec.end());

            if (nestModuleVec.size() > 0) {
                // We found a recently nested element!
                auto recentVec = nestModuleVec.at(0)->getRtVec();
                int graphGroup = rt->getGraphGroup();
                auto graph = mp.getFamilyGraph(rt);
            

                // Removing duplicates
                for (auto& r : *recentVec) {
                    r->push(graphRegion.first);
                    bool foundR = false;
                    for (auto &elePtr : graph->getValueVec()) {
                        if (r->hasRightLTR() && (r->getLeftLTR()->calcOverlap(*elePtr) > 0 || r->getRightLTR()->calcOverlap(*elePtr) > 0)) {
                            foundR = true;
                            break;
                        }
                    }

                    

                    if (!foundR) {
                        delete r;
                        r = nullptr;
                    }
                    else {
                        // delete r;
                        // r = nullptr;
                        r->setGraphGroup(graphGroup);
                        recentNestVec.push_back(r);
                    }
                }

            }
        }
    }

    rtVecPtr->insert(rtVecPtr->end(), recentNestVec.begin(), recentNestVec.end());

}

bool DeepNesting::shouldSkipRegion(int removeStart, int removeLength, int level) {
    // Check recursion depth limit
    if (level >= MAX_NESTING_DEPTH) {
        return true;  // Silently skip without logging
    }
    
    std::pair<int, int> region = std::make_pair(removeStart, removeLength);
    bool shouldSkip = false;
    
    #pragma omp critical
    {
        // Immediately check if we've seen this exact region before
        if (visitedRegions.find(region) != visitedRegions.end()) {
            shouldSkip = true;
        }
        else {
            // Mark as visited BEFORE incrementing count to prevent race conditions
            visitedRegions.insert(region);
            
            // Check visit count for similar regions (with tolerance)
            int similarCount = 0;
            for (const auto& visited : visitedRegions) {
                // Check if regions overlap significantly (>80% overlap)
                int visitedStart = visited.first;
                int visitedLen = visited.second;
                int overlapStart = std::max(removeStart, visitedStart);
                int overlapEnd = std::min(removeStart + removeLength, visitedStart + visitedLen);
                if (overlapEnd > overlapStart) {
                    int overlapLen = overlapEnd - overlapStart;
                    if (overlapLen > removeLength * 0.8 || overlapLen > visitedLen * 0.8) {
                        similarCount++;
                    }
                }
            }
            
            if (similarCount > MAX_SAME_REGION_VISITS) {
                shouldSkip = true;
            }
        }
    }
    
    return shouldSkip;
}

std::vector<ModulePipeline*> DeepNesting::findDeep(std::string &graphSeq, int removeStart, int removeLength, int level) {
    // Global limit check
    bool exceedsLimit = false;
    #pragma omp critical
    {
        totalRegionsProcessed++;
        if (totalRegionsProcessed > MAX_TOTAL_REGIONS) {
            exceedsLimit = true;
        }
    }
    
    if (exceedsLimit) {
        return std::vector<ModulePipeline*>();
    }
    
    // Early termination check
    if (shouldSkipRegion(removeStart, removeLength, level)) {
        return std::vector<ModulePipeline*>();
    }
    
    // Check cache first
    std::pair<int, int> region = std::make_pair(removeStart, removeLength);
    bool foundInCache = false;
    #pragma omp critical
    {
        if (regionCache.find(region) != regionCache.end()) {
            // Found cached result
            foundInCache = true;
        }
    }
    
    if (foundInCache) {
        return std::vector<ModulePipeline*>();
    }
    
    // Optimization: Skip regions that are too large or too small
    if (removeLength > graphSeq.length() / 4 || removeLength < 100) {
        return std::vector<ModulePipeline*>();
    }
    
    std::string cutSeq = graphSeq.substr(0, removeStart) + graphSeq.substr(removeStart + removeLength);
    
    // Early exit if resulting sequence is too small or too similar to original
    if (cutSeq.length() < 1000 || cutSeq.length() > graphSeq.length() * 0.95) {
        return std::vector<ModulePipeline*>();
    }
    
    ModulePipeline *mp = new ModulePipeline{red};
    std::vector<ModulePipeline*> mpVec{mp};

    mp->buildElements(&cutSeq);
    mp->matchElements(icStandard, icRecent, &cutSeq);
    mp->findRTs();
    mp->process(icStandard, &cutSeq, true);

    // Removing complex
    for (auto cpx : *mp->getComplexVec()) {
        delete cpx;
        cpx = nullptr;
    }

    auto rtVec = mp->getRtVec();
    // for (auto &rt : *rtVec) {
    //     delete rt;
    // }
    std::vector<RT*> r;
    
    // Count recently nested elements to limit processing
    int recentNestedCount = 0;
    for (auto &rt : *rtVec) {
        if (rt->getCaseType() == "RecentlyNestedInner") {
            recentNestedCount++;
        }
    }
    
    // Limit the number of recently nested elements processed per level
    const int MAX_RECENT_NESTED_PER_LEVEL = 2;  // Reduced from 5 to 2 for faster processing
    
    int processedCount = 0;
    for (auto &rt : *rtVec) {
        // Searching for the recently nested found inside
        if (rt->getCaseType() == "RecentlyNestedInner") {
            if (processedCount >= MAX_RECENT_NESTED_PER_LEVEL) {
                break;  // Stop processing after limit
            }
            
            // Silently process nests without logging
            
            auto nestModulePipeline = findDeep(cutSeq, rt->getStart(), rt->getSize(), level + 1);
            
            // Check if recursion returned valid results
            if (!nestModulePipeline.empty()) {
                auto nestVec = nestModulePipeline.front()->getRtVec();
                r.insert(r.end(), nestVec->begin(), nestVec->end());
                mpVec.insert(mpVec.end(), nestModulePipeline.begin(), nestModulePipeline.end());
            }
            
            processedCount++;
        }
    }
    rtVec->insert(rtVec->end(), r.begin(), r.end());
    for (auto &rt : *rtVec) {
        rt->expand(removeStart, removeLength);
    }

    rtVec->shrink_to_fit();

    // Cache the results for this region
    #pragma omp critical
    {
        regionCache[region] = *rtVec;
    }

    LtrUtility::removeNests(*rtVec);
    return mpVec;
}