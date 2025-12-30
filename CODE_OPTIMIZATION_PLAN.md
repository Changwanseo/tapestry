# Tapestry Code-Level Optimization Analysis

## Time Breakdown (23m 26s total for 43.7 Mbp assembly)

```
┌─────────────────────────────────┬────────┬──────┬──────────────────┐
│ Step                            │ Time   │ %    │ Parallelized?    │
├─────────────────────────────────┼────────┼──────┼──────────────────┤
│ Read sampling                   │ 6s     │ 0.4% │ ✅ Yes (I/O)     │
│ Read alignment (minimap2)       │ 5m23s  │ 23%  │ ✅ Yes (-t flag) │
│ Contig alignment (minimap2)     │ 1m20s  │ 6%   │ ✅ Yes (-t flag) │
│ Database building               │ 5s     │ 0.4% │ ❌ No            │
│ **Finding neighbours**          │ 13m48s │ 59%  │ ❌❌ **NO!!**     │
│ Processing contigs              │ 2m43s  │ 12%  │ ✅ Yes (Pool)    │
│ Ploidy estimates                │ 1s     │ 0.1% │ ✅ Yes (Pool)    │
│ Report generation               │ 1s     │ 0.1% │ ❌ No            │
└─────────────────────────────────┴────────┴──────┴──────────────────┘
```

## 🔴 CRITICAL BOTTLENECK: `find_neighbours()` Function

**Location**: `tapestry/alignments.py:336-393`

**Current Implementation**:
```python
def find_neighbours(self):
    log.info("Finding neighbouring alignments")
    results = self.get_multi_alignments()  # Fetch ~30k alignments
    
    # 🔴 SINGLE-THREADED LOOP - 37 alignments/sec
    for i, row in tqdm(enumerate(results), total=len(results)):
        update = {...}  # Build update dict
        
        # Sequential array access
        if i > 0 and row[1] == results[i-1][1]:
            update['a_pre_contig'] = results[i-1][2]
            ...
        
        # Batch updates every 1000 records
        if alncount == 1000:
            conn.execute(update_stmt, updates)
            updates = []
```

**Performance Issues**:
1. ❌ Single-threaded iteration (~37 alignments/sec)
2. ❌ List indexing `results[i-1]` is slow for large datasets
3. ❌ Multiple database transactions (31 transactions for 31k alignments)
4. ❌ Python loop overhead for simple comparisons

## Optimization Strategies

### Strategy 1: Vectorize with Pandas/NumPy (RECOMMENDED)
**Impact**: 10-20x faster (13m → 40s-1.5m)

**Approach**: Convert to DataFrame and use vectorized operations

```python
import pandas as pd
import numpy as np

def find_neighbours_vectorized(self):
    log.info("Finding neighbouring alignments (optimized)")
    results = self.get_multi_alignments()
    
    # Convert to DataFrame for vectorized operations
    df = pd.DataFrame(results, columns=['id', 'query', 'contig', 
                                        'qstart', 'qend', 'rlen', 'reversed'])
    
    # Group by read and compute neighbors in one pass
    df['prev_query'] = df['query'].shift(1)
    df['next_query'] = df['query'].shift(-1)
    
    # Vectorized neighbor detection
    df['pre_contig'] = np.where(df['query'] == df['prev_query'], 
                                df['contig'].shift(1), None)
    df['pre_distance'] = np.where(df['query'] == df['prev_query'],
                                  df['qstart'] - df['qend'].shift(1),
                                  df['qstart'] - 1)
    
    df['post_contig'] = np.where(df['query'] == df['next_query'],
                                 df['contig'].shift(-1), None)
    df['post_distance'] = np.where(df['query'] == df['next_query'],
                                   df['qstart'].shift(-1) - df['qend'],
                                   df['rlen'] - df['qend'])
    
    # Handle reversed reads
    mask = df['reversed']
    df.loc[mask, ['pre_contig', 'post_contig']] = \
        df.loc[mask, ['post_contig', 'pre_contig']].values
    df.loc[mask, ['pre_distance', 'post_distance']] = \
        df.loc[mask, ['post_distance', 'pre_distance']].values
    
    # Single bulk update
    updates = df[['id', 'pre_contig', 'pre_distance', 
                  'post_contig', 'post_distance']].to_dict('records')
    
    with self.engine.begin() as conn:
        conn.execute(update_stmt, updates)
```

**Advantages**:
- ✅ 10-20x faster (vectorized NumPy operations)
- ✅ Single database transaction
- ✅ Minimal code changes
- ✅ No additional dependencies (pandas/numpy already indirect deps)

### Strategy 2: Multiprocessing Chunked Processing
**Impact**: 3-4x faster with 4 cores (13m → 3-4m)

**Approach**: Split alignments into chunks, process in parallel

```python
from multiprocessing import Pool
from functools import partial

def process_alignment_chunk(chunk, all_results):
    """Process a chunk of alignments"""
    updates = []
    for local_i, global_i in enumerate(chunk):
        row = all_results[global_i]
        update = {
            'a_id': row[0],
            'a_pre_contig': None, 'a_pre_distance': None,
            'a_post_contig': None, 'a_post_distance': None
        }
        
        # Check previous (need to check chunk boundaries)
        if global_i > 0 and row[1] == all_results[global_i-1][1]:
            update['a_pre_contig'] = all_results[global_i-1][2]
            update['a_pre_distance'] = row[3] - all_results[global_i-1][4]
        else:
            update['a_pre_distance'] = row[3] - 1
        
        # Check next
        if global_i < len(all_results)-1 and row[1] == all_results[global_i+1][1]:
            update['a_post_contig'] = all_results[global_i+1][2]
            update['a_post_distance'] = all_results[global_i+1][3] - row[4]
        else:
            update['a_post_distance'] = row[5] - row[4]
        
        if row[6]:  # Reversed
            pc, pd = update['a_pre_contig'], update['a_pre_distance']
            update['a_pre_contig'], update['a_pre_distance'] = \
                update['a_post_contig'], update['a_post_distance']
            update['a_post_contig'], update['a_post_distance'] = pc, pd
        
        updates.append(update)
    return updates

def find_neighbours_parallel(self, cores=4):
    log.info("Finding neighbouring alignments (parallel)")
    results = self.get_multi_alignments()
    
    # Split into chunks
    chunk_size = len(results) // cores
    chunks = [range(i, min(i+chunk_size, len(results))) 
              for i in range(0, len(results), chunk_size)]
    
    # Process in parallel
    process_func = partial(process_alignment_chunk, all_results=results)
    with Pool(cores) as pool:
        chunk_updates = pool.map(process_func, chunks)
    
    # Flatten and bulk update
    all_updates = [u for chunk in chunk_updates for u in chunk]
    
    update_stmt = ...  # Same as before
    with self.engine.begin() as conn:
        # Batch in chunks of 10000 for large datasets
        for i in range(0, len(all_updates), 10000):
            conn.execute(update_stmt, all_updates[i:i+10000])
```

**Advantages**:
- ✅ Uses existing multicore infrastructure
- ✅ 3-4x speedup
- ✅ Handles chunk boundaries correctly

**Disadvantages**:
- ❌ Still has Python loop overhead
- ❌ More complex than vectorization

### Strategy 3: Database-Level Optimization
**Impact**: 2-3x faster (13m → 4-6m)

**Approach**: Use SQLite window functions and bulk operations

```python
def find_neighbours_sql(self):
    log.info("Finding neighbouring alignments (SQL)")
    
    # Create indexed temp table
    with self.engine.begin() as conn:
        conn.execute(text("""
            CREATE TEMP TABLE temp_neighbors AS
            WITH multi_alns AS (
                SELECT id, query, contig, query_start, query_end, 
                       reversed,
                       LAG(contig) OVER (PARTITION BY query ORDER BY query_start) as pre_contig,
                       LAG(query_end) OVER (PARTITION BY query ORDER BY query_start) as pre_end,
                       LEAD(contig) OVER (PARTITION BY query ORDER BY query_start) as post_contig,
                       LEAD(query_start) OVER (PARTITION BY query ORDER BY query_start) as post_start,
                       (SELECT length FROM reads WHERE name = query) as rlen
                FROM alignments
                WHERE querytype = 'read' 
                  AND (alntype = 'primary' OR alntype = 'supplementary')
                  AND query IN (
                      SELECT query FROM alignments 
                      WHERE querytype='read' 
                        AND (alntype='primary' OR alntype='supplementary')
                      GROUP BY query 
                      HAVING COUNT(*) > 1
                  )
            )
            SELECT 
                id,
                CASE WHEN reversed = 0 THEN pre_contig ELSE post_contig END as pre_contig,
                CASE WHEN reversed = 0 
                     THEN COALESCE(query_start - pre_end, query_start - 1)
                     ELSE COALESCE(post_start - query_end, rlen - query_end) 
                END as pre_distance,
                CASE WHEN reversed = 0 THEN post_contig ELSE pre_contig END as post_contig,
                CASE WHEN reversed = 0 
                     THEN COALESCE(post_start - query_end, rlen - query_end)
                     ELSE COALESCE(query_start - pre_end, query_start - 1)
                END as post_distance
            FROM multi_alns
        """))
        
        # Bulk update from temp table
        conn.execute(text("""
            UPDATE alignments
            SET pre_contig = temp_neighbors.pre_contig,
                pre_distance = temp_neighbors.pre_distance,
                post_contig = temp_neighbors.post_contig,
                post_distance = temp_neighbors.post_distance
            FROM temp_neighbors
            WHERE alignments.id = temp_neighbors.id
        """))
        
        conn.execute(text("DROP TABLE temp_neighbors"))
```

**Advantages**:
- ✅ All processing in SQLite (compiled C code)
- ✅ Uses window functions (LAG/LEAD)
- ✅ Single transaction
- ✅ Leverages database indexes

**Disadvantages**:
- ❌ Complex SQL
- ❌ May hit SQLite limitations on very large datasets

### Strategy 4: Combined Pandas + Cython (MAXIMUM PERFORMANCE)
**Impact**: 20-50x faster (13m → 15-40s)

Use Pandas for vectorization + Cython for hot loops

## Implementation Recommendation

### Phase 1: Quick Win (1-2 hours work)
Implement **Strategy 1 (Pandas vectorization)**

**Steps**:
1. Add pandas import to `tapestry/alignments.py`
2. Replace `find_neighbours()` with vectorized version
3. Test on existing test dataset
4. Verify results match original

**Expected result**: 13m → 40s-1.5m (10-20x faster)

### Phase 2: Additional Speedups (2-3 hours work)
**Database Indexes**:
```python
# Add indexes in table creation (alignments.py:106)
Table('alignments', self.metadata,
    ...existing columns...,
    Index('idx_query_type', 'query', 'querytype'),
    Index('idx_query_alntype', 'query', 'alntype'),
    Index('idx_id', 'id')
)
```

**Expected additional speedup**: 1.5-2x on database queries

### Phase 3: Optimize Other Components (if needed)
1. Parallelize database loading (currently single-threaded)
2. Use memory-mapped I/O for large BAM files
3. Optimize contig processing (already parallelized but could be faster)

## Expected Final Performance

| Configuration | Current | After Phase 1 | After Phase 2 |
|--------------|---------|---------------|---------------|
| 43.7 Mb assembly | 23m 26s | **10-12m** | **8-10m** |
| Time reduction | - | ~60% faster | ~65% faster |

## Code Changes Required

**Minimal changes needed**:
- File: `tapestry/alignments.py`
- Function: `find_neighbours()` (lines 336-393)
- Dependencies: pandas (already indirect dependency via other packages)
- Test coverage: Use existing test_assembly.fasta

## Next Steps

1. ✅ Approve optimization approach
2. Implement Phase 1 (pandas vectorization)
3. Test and benchmark
4. Optionally implement Phase 2 (indexes)
5. Update documentation

