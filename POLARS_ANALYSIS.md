# Polars Performance Analysis

## Test Results

### Original (Loop-based):
- Time: 13m 48s for 31,582 alignments
- Speed: ~38 alignments/second

### Polars v1 (with list comprehensions):
- Time: 13m 54s  
- Speed: ~38 alignments/second
- **NO IMPROVEMENT**

### Polars v2 (optimized DataFrame creation):
- Time: 13m 42s
- Speed: ~38 alignments/second  
- Improvement: Only 12 seconds faster

## Problem Identified

The Polars operations are NOT providing the expected speedup. 

### Why?

The bottleneck is likely NOT in the neighbor-finding logic, but in:

1. **Database UPDATE operations** - Even with batching, updating 31k rows takes time
2. **Data conversion overhead** - Converting between SQL results → DataFrame → dict → SQL

3. **Polars lazy evaluation** - The `.collect()` might be materializing everything at once

## Alternative Approach: NumPy

Since we need raw speed and the operations are simple array manipulations, **NumPy might be better**:

```python
import numpy as np

def find_neighbours_numpy(self):
    results = self.get_multi_alignments()
    
    # Convert to NumPy arrays (FAST!)
    data = np.array(results)
    ids = data[:, 0]
    queries = data[:, 1]
    contigs = data[:, 2]
    qstarts = data[:, 3].astype(int)
    qends = data[:, 4].astype(int)
    rlens = data[:, 5].astype(int)
    reversed_flags = data[:, 6].astype(bool)
    
    # Vectorized neighbor detection (FAST!)
    same_query_prev = queries[:-1] == queries[1:]
    same_query_next = queries[1:] == queries[:-1]
    
    # Pre-allocate arrays
    pre_contigs = np.full(len(queries), None, dtype=object)
    pre_distances = np.zeros(len(queries), dtype=int)
    post_contigs = np.full(len(queries), None, dtype=object)
    post_distances = np.zeros(len(queries), dtype=int)
    
    # Compute distances (vectorized)
    pre_distances[1:] = np.where(same_query_prev, 
                                  qstarts[1:] - qends[:-1],
                                  qstarts[1:] - 1)
    pre_distances[0] = qstarts[0] - 1
    
    post_distances[:-1] = np.where(same_query_next,
                                    qstarts[1:] - qends[:-1],
                                    rlens[:-1] - qends[:-1])
    post_distances[-1] = rlens[-1] - qends[-1]
    
    # Handle reversed (vectorized with np.where)
    ...
```

**Expected speedup:** 10-50x (most operations are pure NumPy)

## Real Issue: Database Updates

Even if we make neighbor-finding instant, the **database UPDATE is the real bottleneck**.

Current: 31,582 individual row updates (even batched)

Better approach:
1. Create temp table with all updates
2. Single UPDATE FROM temp table  
3. Drop temp table

This is what I suggested earlier with the SQL window functions approach!

## Recommendation

**Try SQL window functions approach** - it keeps everything in the database:
- No Python loop overhead
- No data conversion
- Single SQL transaction
- Leverages SQLite's compiled C code

Expected time: **30 seconds to 2 minutes** for the entire operation
