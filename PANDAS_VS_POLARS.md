# Pandas vs Polars for Tapestry Optimization

## Current Status

✅ **Pandas is ALREADY a dependency** (setup.py line 23)
- Used in 2 places: lines 469 and 539 in alignments.py
- Simple DataFrame conversions from SQL results
- Version installed: pandas 2.3.3

## Performance Comparison

### Polars Advantages

| Feature | Pandas | Polars | Winner |
|---------|--------|--------|--------|
| **Speed** | Fast | **2-5x faster** | 🚀 Polars |
| **Memory** | High | **Lower (50-80%)** | 🚀 Polars |
| **Parallelization** | Limited (GIL) | **Native multi-threading** | 🚀 Polars |
| **API** | Implicit | **Explicit & safer** | 🚀 Polars |
| **Lazy evaluation** | No | **Yes (query optimization)** | 🚀 Polars |
| **Installation** | Already installed | Need to add | ⚠️ Pandas |
| **Maturity** | Very mature | Newer (but stable) | ⚠️ Pandas |
| **Learning curve** | Familiar | Different syntax | ⚠️ Pandas |

### Benchmark: Expected Performance for find_neighbours()

**Current**: 13m 48s (single-threaded loop, 31,582 alignments)

| Approach | Time | Speedup | Notes |
|----------|------|---------|-------|
| Current (loop) | 13m 48s | 1x | Baseline |
| Pandas vectorized | 40s - 1.5m | **10-20x** | Good |
| **Polars vectorized** | **15-30s** | **25-50x** | **Best!** |
| Polars + lazy | 10-20s | 40-80x | Maximum |

## Code Comparison

### Option 1: Pandas Implementation

```python
import pandas as pd
import numpy as np

def find_neighbours_pandas(self):
    log.info("Finding neighbouring alignments (pandas)")
    results = self.get_multi_alignments()
    
    # Convert to DataFrame
    df = pd.DataFrame(results, columns=['id', 'query', 'contig', 
                                        'qstart', 'qend', 'rlen', 'reversed'])
    
    # Compute neighbors with shift operations
    df['prev_query'] = df['query'].shift(1)
    df['next_query'] = df['query'].shift(-1)
    
    # Vectorized operations
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
    
    # Bulk update
    updates = df[['id', 'pre_contig', 'pre_distance', 
                  'post_contig', 'post_distance']].to_dict('records')
    
    update_stmt = (self.alignments.update()
                        .where(self.alignments.c.id == bindparam('id'))
                        .values(pre_contig=bindparam('pre_contig'),
                                pre_distance=bindparam('pre_distance'),
                                post_contig=bindparam('post_contig'),
                                post_distance=bindparam('post_distance'))
                   )
    
    with self.engine.begin() as conn:
        conn.execute(update_stmt, updates)
```

**Pros**:
- ✅ Already installed (no new dependency)
- ✅ Familiar syntax
- ✅ Works well with existing code
- ✅ 10-20x speedup

**Cons**:
- ❌ Slower than Polars
- ❌ Higher memory usage
- ❌ Limited parallelization

### Option 2: Polars Implementation (RECOMMENDED)

```python
import polars as pl

def find_neighbours_polars(self):
    log.info("Finding neighbouring alignments (polars)")
    results = self.get_multi_alignments()
    
    # Convert to Polars DataFrame (faster than pandas)
    df = pl.DataFrame({
        'id': [r[0] for r in results],
        'query': [r[1] for r in results],
        'contig': [r[2] for r in results],
        'qstart': [r[3] for r in results],
        'qend': [r[4] for r in results],
        'rlen': [r[5] for r in results],
        'reversed': [r[6] for r in results]
    })
    
    # Polars lazy evaluation for query optimization
    df = df.lazy()
    
    # Compute neighbors with shift (Polars syntax)
    df = df.with_columns([
        pl.col('query').shift(1).alias('prev_query'),
        pl.col('query').shift(-1).alias('next_query'),
        pl.col('contig').shift(1).alias('prev_contig_temp'),
        pl.col('contig').shift(-1).alias('next_contig_temp'),
        pl.col('qend').shift(1).alias('prev_qend'),
        pl.col('qstart').shift(-1).alias('next_qstart')
    ])
    
    # Vectorized neighbor detection (Polars when-then-otherwise)
    df = df.with_columns([
        pl.when(pl.col('query') == pl.col('prev_query'))
          .then(pl.col('prev_contig_temp'))
          .otherwise(None)
          .alias('pre_contig'),
        
        pl.when(pl.col('query') == pl.col('prev_query'))
          .then(pl.col('qstart') - pl.col('prev_qend'))
          .otherwise(pl.col('qstart') - 1)
          .alias('pre_distance'),
        
        pl.when(pl.col('query') == pl.col('next_query'))
          .then(pl.col('next_contig_temp'))
          .otherwise(None)
          .alias('post_contig'),
        
        pl.when(pl.col('query') == pl.col('next_query'))
          .then(pl.col('next_qstart') - pl.col('qend'))
          .otherwise(pl.col('rlen') - pl.col('qend'))
          .alias('post_distance')
    ])
    
    # Handle reversed reads (Polars conditional swapping)
    df = df.with_columns([
        pl.when(pl.col('reversed'))
          .then(pl.col('post_contig'))
          .otherwise(pl.col('pre_contig'))
          .alias('pre_contig_final'),
        
        pl.when(pl.col('reversed'))
          .then(pl.col('post_distance'))
          .otherwise(pl.col('pre_distance'))
          .alias('pre_distance_final'),
        
        pl.when(pl.col('reversed'))
          .then(pl.col('pre_contig'))
          .otherwise(pl.col('post_contig'))
          .alias('post_contig_final'),
        
        pl.when(pl.col('reversed'))
          .then(pl.col('pre_distance'))
          .otherwise(pl.col('post_distance'))
          .alias('post_distance_final')
    ])
    
    # Execute lazy query and convert to dict for SQL update
    df = df.select([
        'id', 
        pl.col('pre_contig_final').alias('pre_contig'),
        pl.col('pre_distance_final').alias('pre_distance'),
        pl.col('post_contig_final').alias('post_contig'),
        pl.col('post_distance_final').alias('post_distance')
    ]).collect()  # Execute lazy query
    
    # Convert to dict for SQL update
    updates = df.to_dicts()
    
    update_stmt = (self.alignments.update()
                        .where(self.alignments.c.id == bindparam('id'))
                        .values(pre_contig=bindparam('pre_contig'),
                                pre_distance=bindparam('pre_distance'),
                                post_contig=bindparam('post_contig'),
                                post_distance=bindparam('post_distance'))
                   )
    
    with self.engine.begin() as conn:
        conn.execute(update_stmt, updates)
```

**Pros**:
- ✅ **2-5x faster than Pandas** (25-50x faster than current)
- ✅ **Lower memory usage** (important for large genomes)
- ✅ **True parallelization** (no GIL limitations)
- ✅ **Lazy evaluation** (query optimization)
- ✅ More explicit API (safer)
- ✅ Modern, actively developed

**Cons**:
- ❌ Adds new dependency (~20 MB)
- ❌ Different syntax (but clearer)
- ❌ Slight learning curve

### Option 3: Hybrid Approach

Keep existing Pandas uses (lines 469, 539) and use Polars ONLY for find_neighbours():

```python
# At top of file
import pandas as pd  # Keep for existing uses
import polars as pl  # Add for optimization

# Use Pandas where already used (no changes)
depths = pd.DataFrame(results)  # Line 469

# Use Polars for bottleneck
def find_neighbours(self):
    # Polars implementation here
```

**Pros**:
- ✅ Best of both worlds
- ✅ Minimal code changes to existing functions
- ✅ Maximum performance where it matters

## Recommendation: Use Polars

### Why Polars?

1. **Performance**: 25-50x speedup vs current (2-5x vs Pandas)
2. **Memory**: Critical for large genome assemblies (>200 Mbp)
3. **Future-proof**: Modern, actively developed
4. **Parallelization**: Native multi-threading (better core utilization)
5. **Small cost**: Just one extra dependency

### Migration Path

**Phase 1**: Replace find_neighbours() with Polars ✅ **DO THIS**
- Biggest impact (59% of runtime)
- Single function change
- Expected: 23m → 8-10m

**Phase 2** (Optional): Keep Pandas for other uses
- Lines 469 and 539 are minimal (not bottlenecks)
- No need to change if working fine

**Phase 3** (Future): Eventually migrate to Polars everywhere
- For consistency
- Better performance across the board

## Implementation

### Step 1: Install Polars

```bash
conda activate tapestry
pip install polars
```

### Step 2: Update setup.py

```python
install_requires=[
    'biopython',
    'intervaltree',
    'jinja2',
    'numpy',
    'pandas',
    'polars',  # ADD THIS
    'plumbum',
    'pysam',
    'sqlalchemy>=1.4.0',
    'tqdm',
]
```

### Step 3: Implement Polars find_neighbours()

Replace function at alignments.py:336-393

### Step 4: Test and Benchmark

```bash
time ./weave -a test_assembly.fasta -r reads.fastq.gz -t TTAGGG -o test_polars -c 4
```

## Expected Results

| Metric | Current | After Polars |
|--------|---------|--------------|
| Runtime | 23m 26s | **8-10m** |
| Memory | Baseline | **60-70%** of baseline |
| CPU usage | 37% (1 core) | **300-400%** (4 cores) |
| Bottleneck | find_neighbours (59%) | Alignment (40%) |

## Conclusion

✅ **YES, switch to Polars for find_neighbours()!**

- Minimal effort (1-2 hours)
- Maximum impact (60-70% faster overall)
- Better memory efficiency
- Future-proof architecture

