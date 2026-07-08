#pragma once

#include "flat_distmx.h"

//////////////////////
// NOTE -- assumes M=L
// Maps 2D coordinates (i, j) to flat 1D distmx
//////////////////////
inline size_t get_k(int i, int j, int L) {
    // Ensure i is the smaller index for upper triangle
    if (i > j) std::swap(i, j); 
    // distmx[M*i + j - i - 1]
    return (size_t)L * i - (i * (i + 1) / 2) + (j - i - 1);
}

class PUUParsers {
public:
    const float DIST_CUTOFF_ANG = 8.0f;
    const uint16_t SID_THRESHOLD;
    const int MIN_DOM_SIZE = 50;

    PUUParsers() : 
        // Convert 8.0A to sid_t: (8*8 * 100) / 16 = 400
        SID_THRESHOLD((uint16_t)((DIST_CUTOFF_ANG * DIST_CUTOFF_ANG * 100) / 16)) {}

    struct Domain {
        int start, end;
    };

    // Recursive function to find domains
void partition(const sid_t* distmx, int L, int start, int end, std::vector<Domain>& results) {
    int currentSize = end - start + 1;

    // Guard: If current segment is too small to be split into two 
    // valid domains, stop recursion and save as a single domain.
    if (currentSize < (2 * MIN_DOM_SIZE)) {
        results.push_back({start, end});
        return;
    }

    int bestSplit = -1;
    float minEnergy = 1e10f;

    // Iterate through all possible split points 's'.
    // The first segment will be [start, s], the second [s+1, end].
    // Both must be >= MIN_DOM_SIZE.
    int loopStart = start + MIN_DOM_SIZE - 1;
    int loopEnd   = end - MIN_DOM_SIZE;

    for (int s = loopStart; s <= loopEnd; ++s) {
        // Validation check for split boundaries
        if (s < 0 || s >= L - 1) continue; 

        float energy = calculateSplitEnergy(distmx, L, start, s, end);
        if (energy < minEnergy) {
            minEnergy = energy;
            bestSplit = s;
        }
    }

    // Threshold check (0.15 is a heuristic; tune as needed)
    if (bestSplit != -1 && minEnergy < 0.15f) {
        partition(distmx, L, start, bestSplit, results);
        partition(distmx, L, bestSplit + 1, end, results);
    } else {
        results.push_back({start, end});
    }
}

private:
    float calculateSplitEnergy(const sid_t* distmx, int L, int start, int split, int end) {
        int interContacts = 0;
        int intraA = 0;
        int intraB = 0;

        for (int i = start; i <= end; ++i) {
            for (int j = i + 1; j <= end; ++j) {
                if (distmx[get_k(i, j, L)] < SID_THRESHOLD) {
                    // Check if contact crosses the split boundary
                    if (i <= split && j > split) {
                        interContacts++;
                    } else if (i <= split && j <= split) {
                        intraA++;
                    } else {
                        intraB++;
                    }
                }
            }
        }

        // Normalized Interaction Energy
        // Avoid division by zero for very small or non-contacting loops
        if (intraA == 0 || intraB == 0) return 1.0f;
        return (float)interContacts / (float)(intraA + intraB);
    }
};