#include "myutils.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

// First draft was by Gemini, then manually edited

using namespace std;

// Structure to hold interval data
struct Interval {
    int start;
    int finish;
    int weight;
    int index; // Original index for tracking (optional)
};

// Comparator for sorting by finish time
bool compareIntervals(const Interval& a, const Interval& b) {
    return a.finish < b.finish;
}

/**
 * Finds the index of the latest non-overlapping interval
 * that finishes before the current interval (i_j).
 * This uses binary search for O(log n) efficiency.
 * * @param sortedIntervals The vector of intervals sorted by finish time.
 * @param j The index of the current interval i_j.
 * @return The index i (0-based) such that i < j and sortedIntervals[i].finish <= sortedIntervals[j].start.
 * Returns -1 if no such non-overlapping interval exists.
 */
int findPredecessor(const vector<Interval>& sortedIntervals, int j) {
    int low = 0;
    int high = j - 1;
    int latest_compatible = -1;
    
    // Binary Search to find the rightmost compatible interval
    while (low <= high) {
        int mid = low + (high - low) / 2;
        
        // Check for compatibility: i_mid finishes before i_j starts
        if (sortedIntervals[mid].finish <= sortedIntervals[j].start) {
            latest_compatible = mid; // Found a compatible one, try for a later one
            low = mid + 1;
        } else {
            high = mid - 1; // Overlaps, need to look earlier
        }
    }
    return latest_compatible;
}

/**
 * Solves the Weighted Interval Scheduling Problem using Dynamic Programming.
 * Time Complexity: O(n log n) due to sorting and n binary searches.
 * Space Complexity: O(n) for DP array and predecessors array.
 * * @param intervals The initial vector of Interval structs.
 * @return The maximum total weight of non-overlapping intervals.
 */
int weightedIntervalScheduling(vector<Interval>& intervals) {
    size_t n = intervals.size();
    if (n == 0) return 0;

    // Step 1: Sort by finish time
    sort(intervals.begin(), intervals.end(), compareIntervals);

    // Step 2: Precompute predecessors
    // predecessor[j] will store the 0-based index of p(j)
    vector<int> predecessor(n);
    for (int j = 0; j < n; ++j) {
        predecessor[j] = findPredecessor(intervals, j);
    }
    
    // Step 3: Dynamic Programming
    // M[j] stores the maximum weight considering intervals 0 to j-1
    // M array size n+1 to handle M[0] base case easily.
    vector<int> M(n + 1);
    M[0] = 0; // Base case: Max weight from 0 intervals is 0

    for (int j = 1; j <= n; ++j) {
        // Interval index in the 0-based sorted array is j-1
        int current_idx = j - 1;
        
        // Case 1: DO NOT include interval i_{j-1}
        int option1 = M[j - 1]; 

        // Case 2: INCLUDE interval i_{j-1}
        int pred_idx = predecessor[current_idx]; // Get the precomputed predecessor index
        
        // If pred_idx is -1, it means no compatible interval exists before it,
        // so M[p(j)] is M[0], which is 0.
        int compatible_weight = (pred_idx != -1) ? M[pred_idx + 1] : 0; 
        
        int option2 = intervals[current_idx].weight + compatible_weight;

        // Take the maximum of the two options
        M[j] = max(option1, option2);
    }

    // M[n] holds the maximum weight considering all n intervals
    return M[n];
}

// --- Function to reconstruct the actual set of intervals (Optional) ---
/**
 * Function to reconstruct the set of selected intervals.
 * * @param M The DP array M[j].
 * @param predecessor The precomputed predecessor array.
 * @param intervals The sorted interval vector.
 * @param j Current index (1-based) in the DP array.
 * @param result The vector to store the selected intervals.
 */
void findSolution(const vector<int>& M, const vector<int>& predecessor, 
                  const vector<Interval>& intervals, int j, vector<Interval>& result) {
    if (j == 0) {
        return;
    }

    // Interval index in the 0-based sorted array is j-1
    int current_idx = j - 1;
    
    // Check if M[j] came from including interval i_{j-1} (Option 2)
    // The compatible weight M[p(j)] is M[predecessor[current_idx] + 1]
    int pred_idx_plus_one = predecessor[current_idx] + 1;
    int compatible_weight = (predecessor[current_idx] != -1) ? M[pred_idx_plus_one] : 0; 
    
    // If the maximum value M[j] is equal to (w_j + M[p(j)]), 
    // it means we included interval i_{j-1}
    if (M[j] == intervals[current_idx].weight + compatible_weight) {
        // Include i_{j-1} in the solution
        result.push_back(intervals[current_idx]);
        // Continue tracking back from its predecessor p(j)
        findSolution(M, predecessor, intervals, pred_idx_plus_one, result);
    } else {
        // Did NOT include i_{j-1}, continue tracking back from M[j-1]
        findSolution(M, predecessor, intervals, j - 1, result);
    }
}


// Example Usage
int example_main() {
    // Intervals: {start, finish, weight}
    vector<Interval> intervals = {
        {1, 4, 20},  // i1
        {3, 5, 20},  // i2
        {0, 6, 10},  // i3
        {5, 7, 70},  // i4
        {3, 8, 30},  // i5
        {5, 9, 40},  // i6
        {6, 10, 60}, // i7
        {8, 11, 50}  // i8
    };

    // Calculate the maximum weight
    int max_weight = weightedIntervalScheduling(intervals);

    cout << "Maximum total weight of non-overlapping intervals: " << max_weight << endl;

    // --- Optional: Reconstruct the actual set of intervals ---
    // Note: The DP array M and predecessor array are not available outside the function scope,
    // a production implementation would return them or structure the code differently.
    // For this demonstration, we'll run the setup steps again to get the data for reconstruction.

    // 1. Sort by finish time (again, since the function modified the vector)
    sort(intervals.begin(), intervals.end(), compareIntervals);

    // 2. Re-compute predecessors
    size_t stn = intervals.size();
    int n = int(stn);
    asserta(size_t(n) == stn);
    vector<int> predecessor(n);
    for (int j = 0; j < n; ++j) {
        predecessor[j] = findPredecessor(intervals, j);
    }

    // 3. Re-compute DP array M (or store it from the main function)
    vector<int> M(n + 1);
    M[0] = 0;
    for (int j = 1; j <= n; ++j) {
        int current_idx = j - 1;
        int option1 = M[j - 1]; 
        int pred_idx_plus_one = predecessor[current_idx] + 1;
        int compatible_weight = (predecessor[current_idx] != -1) ? M[pred_idx_plus_one] : 0; 
        int option2 = intervals[current_idx].weight + compatible_weight;
        M[j] = max(option1, option2);
    }
    
    // 4. Find the solution path
    vector<Interval> solution_set;
    findSolution(M, predecessor, intervals, n, solution_set);

    cout << "\nSelected Optimal Intervals (Start, Finish, Weight):" << endl;
    
    // Output in reverse order since findSolution builds the list backwards
    reverse(solution_set.begin(), solution_set.end()); 
    for (const auto& interval : solution_set) {
        cout << "[" << interval.start << ", " << interval.finish << "], Weight: " << interval.weight << endl;
    }
    
    return 0;
}