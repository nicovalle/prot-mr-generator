#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <utility>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cstdio>

using namespace std;

// LCP function implementations (copied from main program)

bool compute_LCP_naive(const unsigned long n, const unsigned char* T, const unsigned long* SA, unsigned int*& LCP){
    // naive LCP computation (O(n²) time complexity)
    // requires SA previously computed
    
    // a. memory allocation
    LCP = (unsigned int *) malloc((unsigned long)n * sizeof(unsigned int));
    if(LCP == NULL) {
        fprintf(stderr, "Class mers->computeLCP: Cannot allocate memory.\n");
        exit(EXIT_FAILURE);
    }

    // b. LCP computation from SA data
    for (unsigned long iSA = 1; iSA < n; ++iSA) {
        unsigned int l = 0;
        while ((T[SA[iSA]+l]==T[SA[iSA-1]+l]) &&  T[SA[iSA]+l]!='+') ++l;
        LCP[iSA] = l;
    }

    // c. In order to avoid border cases in repeat classification (in other function), we assign LCP[0]=0
    LCP[0]=0;

    return true;
}

bool compute_LCP_kasai(const unsigned long n, const unsigned char* T, const unsigned long* SA, unsigned int*& LCP){
    // Kasai's algorithm for LCP computation (O(n) time complexity)
    // requires SA previously computed
    
    // a. memory allocation
    LCP = (unsigned int *) malloc((unsigned long)n * sizeof(unsigned int));
    if(LCP == NULL) {
        fprintf(stderr, "Class mers->computeLCP_kasai: Cannot allocate memory.\n");
        exit(EXIT_FAILURE);
    }

    // b. Build rank array (inverse of SA): rank[SA[i]] = i
    unsigned long* rank = (unsigned long *) malloc((unsigned long)n * sizeof(unsigned long));
    if(rank == NULL) {
        fprintf(stderr, "Class mers->computeLCP_kasai: Cannot allocate memory for rank array.\n");
        exit(EXIT_FAILURE);
    }
    
    for (unsigned long i = 0; i < n; i++) {
        rank[SA[i]] = i;
    }

    // c. Initialize LCP[0] = 0
    LCP[0] = 0;

    // d. Process suffixes in original text order (not SA order)
    // Key insight: LCP[rank[i]] >= LCP[rank[i-1]] - 1
    unsigned int h = 0;
    for (unsigned long i = 0; i < n; i++) {
        if (rank[i] > 0) {
            unsigned long j = SA[rank[i] - 1];
            // Skip already matched characters using the property above
            while ((i + h < n) && (j + h < n) && 
                   (T[i + h] == T[j + h]) && 
                   (T[i + h] != '+')) {
                h++;
            }
            LCP[rank[i]] = h;
            // Key optimization: LCP[rank[i+1]] >= h-1, so we can start from h-1 next time
            if (h > 0) h--;
        }
    }

    free(rank);
    return true;
}

// Helper function to compare two LCP arrays
bool compareLCPArrays(unsigned int* lcp1, unsigned int* lcp2, unsigned long n) {
    for (unsigned long i = 0; i < n; i++) {
        if (lcp1[i] != lcp2[i]) {
            return false;
        }
    }
    return true;
}

// Helper function to compute checksum of LCP array
unsigned long computeLCPChecksum(unsigned int* lcp, unsigned long n) {
    unsigned long checksum = 0;
    for (unsigned long i = 0; i < n; i++) {
        checksum += lcp[i];
        checksum = (checksum << 1) | (checksum >> 31); // Simple hash
    }
    return checksum;
}

// Helper function to print LCP array differences
void printLCPDiff(unsigned int* lcp1, unsigned int* lcp2, unsigned long n, const unsigned long* SA, const unsigned char* T) {
    cout << "LCP Array Differences:" << endl;
    bool foundDiff = false;
    for (unsigned long i = 0; i < n; i++) {
        if (lcp1[i] != lcp2[i]) {
            foundDiff = true;
            cout << "  Index " << i << ": naive=" << lcp1[i] << ", kasai=" << lcp2[i];
            if (i > 0 && SA != NULL && T != NULL) {
                cout << " (SA[" << i-1 << "]=" << SA[i-1] << ", SA[" << i << "]=" << SA[i] << ")";
            }
            cout << endl;
        }
    }
    if (!foundDiff) {
        cout << "  No differences found." << endl;
    }
}

// Simple SA computation for testing (naive approach)
void computeSA_naive(const unsigned char* T, unsigned long n, unsigned long* SA) {
    // Create array of indices
    vector<pair<string, unsigned long>> suffixes;
    for (unsigned long i = 0; i < n; i++) {
        string suffix = "";
        for (unsigned long j = i; j < n; j++) {
            suffix += (char)T[j];
        }
        suffixes.push_back(make_pair(suffix, i));
    }
    
    // Sort by suffix string
    sort(suffixes.begin(), suffixes.end());
    
    // Fill SA array
    for (unsigned long i = 0; i < n; i++) {
        SA[i] = suffixes[i].second;
    }
}

// Test case runner
bool runTestCase(const string& testName, const string& text, const unsigned long* SA, unsigned long n) {
    cout << "\n=== Test: " << testName << " ===" << endl;
    cout << "Text length: " << n << endl;
    
    // Convert string to unsigned char array
    unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
    for (unsigned long i = 0; i < n; i++) {
        T[i] = (unsigned char)text[i];
    }
    
    // Compute LCP with both algorithms
    unsigned int* lcp_naive = NULL;
    unsigned int* lcp_kasai = NULL;
    
    clock_t start, end;
    
    start = clock();
    bool ok1 = compute_LCP_naive(n, T, SA, lcp_naive);
    end = clock();
    double naive_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    start = clock();
    bool ok2 = compute_LCP_kasai(n, T, SA, lcp_kasai);
    end = clock();
    double kasai_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    if (!ok1 || !ok2) {
        cout << "ERROR: Failed to compute LCP arrays" << endl;
        free(T);
        if (lcp_naive) free(lcp_naive);
        if (lcp_kasai) free(lcp_kasai);
        return false;
    }
    
    // Compare results
    bool match = compareLCPArrays(lcp_naive, lcp_kasai, n);
    unsigned long checksum_naive = computeLCPChecksum(lcp_naive, n);
    unsigned long checksum_kasai = computeLCPChecksum(lcp_kasai, n);
    
    cout << "Naive time: " << naive_time << "s" << endl;
    cout << "Kasai time: " << kasai_time << "s" << endl;
    if (naive_time > 0) {
        cout << "Speedup: " << (naive_time / kasai_time) << "x" << endl;
    }
    cout << "Checksum (naive): " << checksum_naive << endl;
    cout << "Checksum (kasai): " << checksum_kasai << endl;
    
    if (!match) {
        cout << "FAILED: LCP arrays do not match!" << endl;
        printLCPDiff(lcp_naive, lcp_kasai, n, SA, T);
    } else {
        cout << "PASSED: LCP arrays match perfectly" << endl;
    }
    
    free(T);
    free(lcp_naive);
    free(lcp_kasai);
    
    return match;
}

int main() {
    cout << "========================================" << endl;
    cout << "LCP Algorithm Comparison Tests" << endl;
    cout << "Standalone Test Executable" << endl;
    cout << "========================================" << endl;
    
    int passed = 0;
    int total = 0;
    
    // Test Case 1: Simple string without separators
    {
        string text = "banana$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("Simple string (banana)", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Test Case 2: String with separators (protein-like)
    {
        string text = "ABCD+EFAB+GHAB+$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("String with separators", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Test Case 3: Repeated patterns
    {
        string text = "AAAA+BBBB+AAAA+$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("Repeated patterns", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Test Case 4: Single character
    {
        string text = "A$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("Single character", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Test Case 5: All same characters
    {
        string text = "AAAAA$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("All same characters", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Test Case 6: Protein-like sequence
    {
        string text = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("Protein-like sequence", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Test Case 7: Edge case - minimal string
    {
        string text = "$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("Minimal string (only terminator)", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Test Case 8: Complex pattern with multiple repeats
    {
        string text = "ABABAB+CDCDCD+ABABAB+$";
        unsigned long n = text.length();
        unsigned long* SA = (unsigned long*)malloc(n * sizeof(unsigned long));
        unsigned char* T = (unsigned char*)malloc(n * sizeof(unsigned char));
        for (unsigned long i = 0; i < n; i++) {
            T[i] = (unsigned char)text[i];
        }
        computeSA_naive(T, n, SA);
        total++;
        if (runTestCase("Complex pattern with multiple repeats", text, SA, n)) {
            passed++;
        }
        free(SA);
        free(T);
    }
    
    // Summary
    cout << "\n========================================" << endl;
    cout << "Test Summary" << endl;
    cout << "========================================" << endl;
    cout << "Passed: " << passed << "/" << total << endl;
    if (passed == total) {
        cout << "All tests PASSED!" << endl;
        return 0;
    } else {
        cout << "Some tests FAILED!" << endl;
        return 1;
    }
}
