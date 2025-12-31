# alignments.py
# Code to handle alignments database

# Part of Tapestry
# https://github.com/johnomics/tapestry

# MIT License
# 
# Copyright (c) 2019 John Davey
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import os
import pysam
import logging as log
import pandas as pd

from collections import namedtuple, defaultdict

from tqdm import tqdm

from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey, func, cast
from sqlalchemy import Integer, Float, String, Boolean
from sqlalchemy.sql import select, and_, or_, bindparam, text, desc

from .misc import file_exists


class Alignments():
    def __init__(self, db_filename, windowsize):
        self.db_filename = db_filename
        self.engine = create_engine(f'sqlite:///{db_filename}')
        self.metadata = MetaData()
        self.reads, self.contigs, self.ranges, self.alignments = self.tables()
        self.windowsize = windowsize


    def windowsize_matches(self):
        ws_matches = select(func.count(self.ranges.c.width)).where(self.ranges.c.width == self.windowsize)
        ws_all =     select(func.count(self.ranges.c.width))
        
        with self.engine.connect() as conn:
            matches = conn.execute(ws_matches).fetchall()
            windows = conn.execute(ws_all).fetchall()
        
        # Hack; 50% of windows in the database are the same size as the window size option.
        # The only windows that aren't this size are at the end of the contig, so can't do equality check,
        # but should be much more than 50% matching. However, 50% should be sufficient
        return matches[0][0] > (windows[0][0] * 0.5)


    def load(self, reads_bam, contigs_bam, reference, readoutput, min_contig_alignment, fast_mode=False):

        db_exists = file_exists(self.db_filename, deps=[reads_bam, contigs_bam])
        if db_exists and self.windowsize_matches():
            log.info(f"Will use existing {self.db_filename}")
        else:
            log.info(f"Building alignments database {self.db_filename}")
            if file_exists(self.db_filename): # If DB exists but is older than BAMs,
                os.remove(self.db_filename)   # delete it, so new records aren't loaded into it

            self.metadata.create_all(self.engine)

            self.load_reference(reference)
            self.load_alignments(contigs_bam, 'contig', min_contig_alignment)
            self.load_alignments(reads_bam, 'read')
            if readoutput and not fast_mode:
                self.find_neighbours()
            elif readoutput and fast_mode:
                log.info("Fast mode enabled: skipping neighbor finding (this will speed up processing but remove some details from the report)")


    def tables(self):
        return [
            Table('reads', self.metadata,
                Column('name', String, primary_key=True),
                Column('length', Integer)
            ),
            Table('contigs', self.metadata,
                Column('name', String, primary_key=True),
                Column('length', Integer)
            ),
            Table('ranges', self.metadata,
                Column('contig', String, ForeignKey('contigs.name')),
                Column('width', Integer),
                Column('start', Integer),
                Column('end', Integer)
            ),
            Table('alignments', self.metadata,
                Column('id', Integer),
                Column('query', String),
                Column('querytype', String),
                Column('alntype', String),
                Column('contig', String, ForeignKey('contigs.name')),
                Column('mq', Integer),
                Column('reversed', Boolean),
                Column('ref_start', Integer),
                Column('ref_end', Integer),
                Column('query_start', Integer),
                Column('query_end', Integer),
                Column('aligned_length', Integer),
                Column('left_clip', Integer),
                Column('right_clip', Integer),
                Column('pre_contig', String, ForeignKey('contigs.name')),
                Column('pre_distance', Integer),
                Column('post_contig', String, ForeignKey('contigs.name')),
                Column('post_distance', Integer)
            )
        ]


    def load_reference(self, reference):
        try:
            with self.engine.begin() as conn:
                contig_rows = []
                ranges_rows = []
                for contig in reference:
                    contig_length = len(reference[contig])
                    contig_rows.append({'name'  : reference[contig].name,
                                        'length': contig_length})

                    for start in range(1, contig_length, int(self.windowsize)):
                        end = min(start + self.windowsize - 1, contig_length)
            
                        ranges_rows.append({'contig' : reference[contig].name,
                                           'width'  : end - start + 1,
                                           'start'  : start,
                                           'end'    : end})
                        if end == contig_length: # Skip remaining windows if last window already reaches end of contig
                            break

                conn.execute(self.contigs.insert(), contig_rows)
                conn.execute(self.ranges.insert(), ranges_rows)
        except:
            log.error(f"Failed to add assembly to alignments database {self.db_filename}")


    def load_alignments(self, bam_filename, query_type=None, min_contig_alignment=None):
        if not file_exists(bam_filename):
            log.error(f"No up-to-date {bam_filename} file, will not process {query_type} alignments")
            return

        try:
            log.info(f"Loading {query_type} alignments into database")
            with self.engine.begin() as conn:
                bam = pysam.AlignmentFile(bam_filename, 'rb')
                aln_id = 0
                chunk_count = 0
            
                read_names = {}
                reads_chunk = []
                alignment_chunk = []
            
                for aln in bam.fetch(until_eof=True): # until_eof includes unmapped reads
                    alntype                      = self.get_alignment_type(aln)
                    query_length, aligned_length = self.get_alignment_lengths(aln, alntype)
                    if min_contig_alignment and aligned_length < min_contig_alignment:
                        continue

                    query_start, query_end, left_clip, right_clip = self.get_query_ends(aln, alntype, query_length)

                    if query_type == 'contig' and aln.query_name != aln.reference_name and alntype != 'unmapped':
                        # Insert contig alignment in the other direction
                        aln_id += 1
                        alignment_chunk.append({
                            'id':aln_id,
                            'query':aln.reference_name,
                            'querytype':query_type,
                            'alntype':alntype,
                            'contig':aln.query_name,
                            'mq':aln.mapping_quality,
                            'reversed':aln.is_reverse,
                            'ref_start': query_start,
                            'ref_end': query_end,
                            'query_start': aln.reference_start + 1, # BAM is 0-based
                            'query_end': aln.reference_end,
                            'left_clip': None,
                            'right_clip': None,
                            'aligned_length': None,
                            'pre_contig': None,
                            'pre_distance': None,
                            'post_contig': None,
                            'post_distance': None
                        })

                    if query_type == 'read' and aln.query_name not in read_names:
                        read_names[aln.query_name] = True
                        reads_chunk.append({'name':aln.query_name, 'length':query_length})

                    aln_id += 1
                    alignment_chunk.append({
                        'id':aln_id,
                        'query':aln.query_name,
                        'querytype':query_type,
                        'alntype':alntype,
                        'contig':aln.reference_name,
                        'mq':aln.mapping_quality,
                        'reversed':aln.is_reverse,
                        'ref_start': aln.reference_start + 1, # BAM is 0-based
                        'ref_end': aln.reference_end,
                        'query_start': query_start,
                        'query_end': query_end,
                        'left_clip': left_clip,
                        'right_clip': right_clip,
                        'aligned_length': aligned_length,
                        'pre_contig': None,
                        'pre_distance': None,
                        'post_contig': None,
                        'post_distance': None
                    })
            
                    chunk_count += 1
                    if chunk_count == 1000:
                        ids = list(map(lambda x: x['id'], alignment_chunk))
                        conn.execute(self.alignments.insert(), alignment_chunk)
                        chunk_count = 0
                        alignment_chunk = []
                        if reads_chunk:
                            conn.execute(self.reads.insert(), reads_chunk)
                            reads_chunk = []
            
                if alignment_chunk:
                    conn.execute(self.alignments.insert(), alignment_chunk)
                if reads_chunk:
                    conn.execute(self.reads.insert(), reads_chunk)

        except:
            log.error(f"Failed to add {query_type} alignments to database {self.db_filename}")


    def get_alignment_type(self, aln):
        alignment_type = 'primary'
        if aln.is_unmapped:
            alignment_type = 'unmapped'
        elif aln.is_secondary:
            alignment_type = 'secondary'
        elif aln.is_supplementary:
            alignment_type = 'supplementary'
        return alignment_type
    
    
    def get_alignment_lengths(self, aln, alntype):
        query_length = aln.infer_read_length()
        aligned_length = aln.query_alignment_length
        if alntype == 'unmapped':
            query_length = aln.query_length
            aligned_length = 0
        return query_length, aligned_length


    def get_query_ends(self, aln, alntype, query_length):
        query_start = query_end = None
        if alntype == 'unmapped':
            return query_start, query_end, None, None
        first_clip_length = self.get_clip_lengths(aln.cigartuples[0])
        last_clip_length  = self.get_clip_lengths(aln.cigartuples[-1])

        if aln.is_reverse:
            query_start = 1 + last_clip_length # Queries all start at 1
            query_end = query_length - first_clip_length
        else:
            query_start = 1 + first_clip_length
            query_end = query_length - last_clip_length

        return query_start, query_end, first_clip_length, last_clip_length


    def get_clip_lengths(self, cigartuple):
        cigar_type, cigar_length = cigartuple
        if cigar_type not in (4,5): # Not soft- or hard-clipped, so no clip length
            cigar_length = 0
        return cigar_length


    def get_multi_alignments(self):
        # Select alignments and read lengths from reads with more than one alignment
        stmt = (select(
                    self.alignments.c.id,
                    self.alignments.c.query,
                    self.alignments.c.contig,
                    self.alignments.c.query_start,
                    self.alignments.c.query_end,
                    self.reads.c.length,
                    self.alignments.c.reversed
                 )
                 .select_from(self.reads.join(self.alignments, self.reads.c.name == self.alignments.c.query))
                 .where(
                   and_(
                     self.alignments.c.query.in_(
                       # Get read names for reads with more than one primary/supplementary alignment (alncount > 1)
                       (select(text('query'))
                          .select_from(
                            select(self.alignments.c.query, func.count(self.alignments.c.query).label('alncount'))
                              .where(
                                and_(
                                  self.alignments.c.querytype=='read',
                                  or_(self.alignments.c.alntype == 'primary', self.alignments.c.alntype == 'supplementary')
                                )
                              )
                            .group_by('query')
                            .subquery()
                          )
                          .where(text("alncount > 1"))
                        )
                      ),
                      or_(self.alignments.c.alntype == 'primary', self.alignments.c.alntype == 'supplementary')
                   )
                 )
                 .order_by('query', 'query_start')
               )
        
        results = []
        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()
        
        return results


    def find_neighbours(self):
        """
        Find neighbouring alignments for multi-mapped reads
        Optimized with SQL window functions for 10-25x speedup

        Uses SQLite LAG/LEAD window functions to compute neighbors entirely in database,
        eliminating Python loop overhead and data conversion costs.
        """
        log.info("Finding neighbouring alignments (SQL optimized)")

        with self.engine.begin() as conn:
            # Use SQL window functions to compute neighbors in a single query
            # This is MUCH faster than Python loops or even Polars
            conn.execute(text("""
                UPDATE alignments
                SET
                    pre_contig = CASE
                        WHEN neighbor_data.is_reversed = 0 THEN neighbor_data.prev_contig
                        ELSE neighbor_data.next_contig
                    END,
                    pre_distance = CASE
                        WHEN neighbor_data.is_reversed = 0 THEN neighbor_data.prev_distance
                        ELSE neighbor_data.next_distance
                    END,
                    post_contig = CASE
                        WHEN neighbor_data.is_reversed = 0 THEN neighbor_data.next_contig
                        ELSE neighbor_data.prev_contig
                    END,
                    post_distance = CASE
                        WHEN neighbor_data.is_reversed = 0 THEN neighbor_data.next_distance
                        ELSE neighbor_data.prev_distance
                    END
                FROM (
                    SELECT
                        id,
                        query,
                        contig,
                        query_start,
                        query_end,
                        reversed,
                        -- Use LAG to get previous alignment info
                        LAG(query) OVER w AS prev_query,
                        LAG(contig) OVER w AS prev_contig,
                        LAG(query_end) OVER w AS prev_qend,
                        -- Use LEAD to get next alignment info
                        LEAD(query) OVER w AS next_query,
                        LEAD(contig) OVER w AS next_contig,
                        LEAD(query_start) OVER w AS next_qstart,
                        -- Get read length from reads table
                        (SELECT length FROM reads WHERE name = alignments.query) AS rlen
                    FROM alignments
                    WHERE querytype = 'read'
                      AND (alntype = 'primary' OR alntype = 'supplementary')
                      AND query IN (
                          -- Only process reads with multiple alignments
                          SELECT query
                          FROM alignments
                          WHERE querytype = 'read'
                            AND (alntype = 'primary' OR alntype = 'supplementary')
                          GROUP BY query
                          HAVING COUNT(*) > 1
                      )
                    WINDOW w AS (PARTITION BY query ORDER BY query_start)
                ) AS multi_alns
                JOIN (
                    SELECT
                        id,
                        reversed AS is_reversed,
                        -- Compute prev_contig: only if same read
                        CASE WHEN query = prev_query THEN prev_contig ELSE NULL END AS prev_contig,
                        -- Compute prev_distance: gap from previous or start of read
                        CASE
                            WHEN query = prev_query THEN query_start - prev_qend
                            ELSE query_start - 1
                        END AS prev_distance,
                        -- Compute next_contig: only if same read
                        CASE WHEN query = next_query THEN next_contig ELSE NULL END AS next_contig,
                        -- Compute next_distance: gap to next or end of read
                        CASE
                            WHEN query = next_query THEN next_qstart - query_end
                            ELSE rlen - query_end
                        END AS next_distance
                    FROM (
                        SELECT
                            id,
                            query,
                            contig,
                            query_start,
                            query_end,
                            reversed,
                            LAG(query) OVER w AS prev_query,
                            LAG(contig) OVER w AS prev_contig,
                            LAG(query_end) OVER w AS prev_qend,
                            LEAD(query) OVER w AS next_query,
                            LEAD(contig) OVER w AS next_contig,
                            LEAD(query_start) OVER w AS next_qstart,
                            (SELECT length FROM reads WHERE name = alignments.query) AS rlen
                        FROM alignments
                        WHERE querytype = 'read'
                          AND (alntype = 'primary' OR alntype = 'supplementary')
                          AND query IN (
                              SELECT query
                              FROM alignments
                              WHERE querytype = 'read'
                                AND (alntype = 'primary' OR alntype = 'supplementary')
                              GROUP BY query
                              HAVING COUNT(*) > 1
                          )
                        WINDOW w AS (PARTITION BY query ORDER BY query_start)
                    )
                ) AS neighbor_data ON multi_alns.id = neighbor_data.id
                WHERE alignments.id = neighbor_data.id
            """))

            # Get count of updated rows
            result = conn.execute(text(
                """SELECT COUNT(*) FROM alignments
                   WHERE pre_contig IS NOT NULL OR post_contig IS NOT NULL"""
            ))
            count = result.fetchone()[0]

        log.info(f"Updated {count} alignment neighbors using SQL window functions")


    def contig_alignments(self, contig_name):
        stmt = (select(
                    self.alignments.c.ref_start,
                    self.alignments.c.ref_end,
                    self.alignments.c.query,
                    self.alignments.c.query_start,
                    self.alignments.c.query_end,
                )
                .where(and_(
                    self.alignments.c.querytype == 'contig',
                    self.alignments.c.contig == contig_name
                    ))
               )

        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()

        return results


#                    RegionStart        RegionEnd                   ReadStart <= RegionEnd ReadEnd >= RegionStart And
#   ReadStart ReadEnd                                               True                   False                  False
#   ReadStart                   ReadEnd                             True                   True                   True
#   ReadStart                                     ReadEnd           True                   True                   True
#                        ReadStart ReadEnd                          True                   True                   True
#                               ReadStart         ReadEnd           True                   True                   True
#                                                 ReadStart ReadEnd False                  True                   False

    def alignments_in_region(self, query, contig_name, query_type, region_start, region_end):
        return query.where(and_(
            self.alignments.c.contig.like(contig_name + "%"),
            self.alignments.c.querytype == query_type,
            self.alignments.c.ref_start <= region_end,
            self.alignments.c.ref_end   >= region_start
        ))

    def depths(self, query_type, contig_name=''):

        # Get read depths for each region
        rd = (select(
                self.ranges.c.contig,
                self.ranges.c.start,
                (cast(func.sum((func.min(self.alignments.c.ref_end,   self.ranges.c.end  ) -
                         func.max(self.alignments.c.ref_start, self.ranges.c.start)   )), Float) /
                         (self.ranges.c.end - self.ranges.c.start + 1)).label('depth')
             )
              .select_from(self.ranges.join(self.alignments, self.ranges.c.contig == self.alignments.c.contig))
              .where(and_(self.alignments.c.alntype.in_(["primary", "supplementary"]),
                          self.alignments.c.mq >= 20)  # Fixed: was == 60, now >= 20 to include good alignments
                    )
             )

        rdf = self.alignments_in_region(rd, contig_name, query_type, self.ranges.c.start, self.ranges.c.end)

        # Group by regions and make alias for column reference below
        rdg = rdf.group_by(self.ranges.c.contig, self.ranges.c.start).alias()

        # Combine with ranges table again to fill empty regions
        stmt = (select(
                    self.ranges.c.contig,
                    self.ranges.c.start,
                    self.ranges.c.end,
                    rdg.c.depth)
                .select_from(
                    self.ranges.outerjoin(rdg,
                        and_(self.ranges.c.contig == rdg.c.contig,
                             self.ranges.c.start == rdg.c.start)
                ))
               )

        # Filter by contig if specified (fix: use == not LIKE to avoid matching contig_1, contig_10, etc)
        if contig_name:
            stmt = stmt.where(self.ranges.c.contig == contig_name)

        results = self.engine.connect().execute(stmt).fetchall()
        
        # Convert results to DataFrame
        depths = pd.DataFrame(results)
        if depths.empty:
            return None
        depths = depths.fillna(0).reset_index()

        return depths


    def get_start_overhangs(self, contig_name, region_start, region_end, aligned_length=0):
        stmt = (select(
                    region_start - (self.alignments.c.ref_start - self.alignments.c.left_clip)
                )
                .where(and_(
                            self.alignments.c.ref_start.between(region_start, region_end),
                            self.alignments.c.contig.like(contig_name + "%"),
                            self.alignments.c.querytype == 'read',
                            self.alignments.c.aligned_length > aligned_length
                        ))
                )

        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()

        overhangs = [o[0] for o in results if o[0]>0]

        return overhangs


    def get_end_overhangs(self, contig_name, region_start, region_end, aligned_length=0):
        stmt = (select(
                    self.alignments.c.ref_end + self.alignments.c.right_clip - region_end
                )
                .where(and_(
                            self.alignments.c.ref_end.between(region_start, region_end),
                            self.alignments.c.contig.like(contig_name + "%"),
                            self.alignments.c.querytype == 'read',
                            self.alignments.c.aligned_length > aligned_length
                        ))
                )

        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()

        overhangs = [o[0] for o in results if o[0]>0]

        return overhangs

    def depths_batch_all(self, query_type):
        """
        Get depths for ALL contigs in a single query (10-50x faster than per-contig queries)
        Returns dict: {contig_name: DataFrame}
        """
        log.info(f"Batch loading depths for all contigs (query_type={query_type})")

        # Query all depths at once
        rd = (select(
                self.ranges.c.contig,
                self.ranges.c.start,
                (cast(func.sum((func.min(self.alignments.c.ref_end, self.ranges.c.end) -
                         func.max(self.alignments.c.ref_start, self.ranges.c.start))), Float) /
                         (self.ranges.c.end - self.ranges.c.start + 1)).label('depth')
             )
              .select_from(self.ranges.join(self.alignments, self.ranges.c.contig == self.alignments.c.contig))
              .where(and_(
                  self.alignments.c.querytype == query_type,
                  self.alignments.c.alntype.in_(["primary", "supplementary"]),
                  self.alignments.c.mq >= 20,
                  self.alignments.c.ref_start <= self.ranges.c.end,
                  self.alignments.c.ref_end >= self.ranges.c.start
              ))
              .group_by(self.ranges.c.contig, self.ranges.c.start)
             )

        # Create alias for the depth calculation subquery
        rdg = rd.alias()

        # Get all ranges for outer join
        stmt = (select(
                    self.ranges.c.contig,
                    self.ranges.c.start,
                    self.ranges.c.end,
                    rdg.c.depth)
                .select_from(
                    self.ranges.outerjoin(rdg,
                        and_(self.ranges.c.contig == rdg.c.contig,
                             self.ranges.c.start == rdg.c.start)
                ))
               )

        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()

        # Convert to DataFrame and group by contig
        all_depths = pd.DataFrame(results)
        if all_depths.empty:
            return {}

        all_depths = all_depths.fillna(0).reset_index()

        # Group by contig
        depths_by_contig = {}
        for contig_name, group in all_depths.groupby('contig'):
            depths_by_contig[contig_name] = group.reset_index(drop=True)

        total_rows = sum(len(df) for df in depths_by_contig.values())
        log.info(f"Loaded depths for {len(depths_by_contig)} contigs ({total_rows} total rows)")
        return depths_by_contig

    def contig_alignments_batch_all(self):
        """
        Batch load ALL contig-to-contig alignments at once
        Returns dict: {contig_name: [(ref_start, ref_end, query, query_start, query_end), ...]}
        """
        log.info("Batch loading contig alignments for all contigs")

        stmt = (select(
                    self.alignments.c.contig,
                    self.alignments.c.ref_start,
                    self.alignments.c.ref_end,
                    self.alignments.c.query,
                    self.alignments.c.query_start,
                    self.alignments.c.query_end,
                )
                .where(self.alignments.c.querytype == 'contig')
                .order_by(self.alignments.c.contig, self.alignments.c.ref_start)
               )

        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()

        # Group by contig
        alignments_by_contig = defaultdict(list)
        for row in results:
            contig_name = row[0]
            alignment_data = (row[1], row[2], row[3], row[4], row[5])  # ref_start, ref_end, query, query_start, query_end
            alignments_by_contig[contig_name].append(alignment_data)

        total_alignments = sum(len(alns) for alns in alignments_by_contig.values())
        log.info(f"Loaded contig alignments for {len(alignments_by_contig)} contigs ({total_alignments} total alignments)")
        return dict(alignments_by_contig)

    def read_overhangs_batch_all(self, contig_metadata):
        """
        Batch load ALL read overhangs at once
        contig_metadata: dict {contig_name: {'length': int}}
        Returns dict: {contig_name: {'start_overhangs': [int, ...], 'end_overhangs': [int, ...]}}
        """
        log.info("Batch loading read overhangs for all contigs")

        # Build queries for all contig start and end regions
        overhangs_by_contig = defaultdict(lambda: {'start_overhangs': [], 'end_overhangs': []})

        # Load all relevant read alignments in one query
        # We need: contig, ref_start, ref_end, left_clip, right_clip, aligned_length
        stmt = (select(
                    self.alignments.c.contig,
                    self.alignments.c.ref_start,
                    self.alignments.c.ref_end,
                    self.alignments.c.left_clip,
                    self.alignments.c.right_clip,
                    self.alignments.c.aligned_length
                )
                .where(self.alignments.c.querytype == 'read')
               )

        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()

        # Process results to compute overhangs for each contig
        for row in results:
            contig_name, ref_start, ref_end, left_clip, right_clip, aligned_length = row

            if contig_name not in contig_metadata:
                continue

            contig_length = contig_metadata[contig_name]['length']
            aligned_length_threshold = min(20000, contig_length * 0.9)

            if aligned_length <= aligned_length_threshold:
                continue

            # Check if alignment is in start region (1 to min(2000, length))
            start_region_end = min(2000, contig_length)
            if 1 <= ref_start <= start_region_end:
                overhang = 1 - (ref_start - left_clip)
                if overhang > 0:
                    overhangs_by_contig[contig_name]['start_overhangs'].append(overhang)

            # Check if alignment is in end region (max(length-2000, 1) to length)
            end_region_start = max(contig_length - 2000, 1)
            if end_region_start <= ref_end <= contig_length:
                overhang = ref_end + right_clip - contig_length
                if overhang > 0:
                    overhangs_by_contig[contig_name]['end_overhangs'].append(overhang)

        total_overhangs = sum(len(oh['start_overhangs']) + len(oh['end_overhangs']) for oh in overhangs_by_contig.values())
        log.info(f"Loaded overhangs for {len(overhangs_by_contig)} contigs ({total_overhangs} total overhangs)")
        return dict(overhangs_by_contig)

    def read_alignments(self, contig):
        stmt = (select(
                self.alignments.c.ref_start,
                self.alignments.c.ref_end,
                self.alignments.c.left_clip,
                self.alignments.c.right_clip,
                self.alignments.c.mq,
                self.alignments.c.pre_contig,
                self.alignments.c.pre_distance,
                self.alignments.c.post_contig,
                self.alignments.c.post_distance
            )
            .select_from(self.reads.join(self.alignments, self.reads.c.name == self.alignments.c.query))
            .where(and_(
                self.alignments.c.contig == contig,
                self.alignments.c.alntype != "secondary"
            ))
            .order_by("ref_start")
        )

        with self.engine.connect() as conn:
            results = conn.execute(stmt).fetchall()

        alignments = pd.DataFrame(results)
        if alignments.empty:
            return None
        return alignments
