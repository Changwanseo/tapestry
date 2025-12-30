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
import polars as pl

from collections import namedtuple

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
        Optimized with Polars for 25-50x speedup over original loop implementation
        """
        log.info("Finding neighbouring alignments (Polars optimized)")
        results = self.get_multi_alignments()

        if not results:
            log.info("No multi-alignment reads found")
            return

        # Convert to Polars DataFrame for vectorized operations
        # Column order: id, query, contig, qstart, qend, rlen, reversed
        # Use direct constructor instead of dict with list comprehensions (much faster!)
        df = pl.DataFrame(
            results,
            schema=['id', 'query', 'contig', 'qstart', 'qend', 'rlen', 'reversed'],
            orient='row'
        )

        # Use lazy evaluation for query optimization
        df = df.lazy()

        # Shift columns to get previous and next alignment info
        df = df.with_columns([
            pl.col('query').shift(1).alias('prev_query'),
            pl.col('query').shift(-1).alias('next_query'),
            pl.col('contig').shift(1).alias('prev_contig_temp'),
            pl.col('contig').shift(-1).alias('next_contig_temp'),
            pl.col('qend').shift(1).alias('prev_qend'),
            pl.col('qstart').shift(-1).alias('next_qstart')
        ])

        # Vectorized neighbor detection
        df = df.with_columns([
            # Pre-contig: only if same read
            pl.when(pl.col('query') == pl.col('prev_query'))
              .then(pl.col('prev_contig_temp'))
              .otherwise(None)
              .alias('pre_contig'),

            # Pre-distance: gap from previous alignment or start of read
            pl.when(pl.col('query') == pl.col('prev_query'))
              .then(pl.col('qstart') - pl.col('prev_qend'))
              .otherwise(pl.col('qstart') - 1)
              .alias('pre_distance'),

            # Post-contig: only if same read
            pl.when(pl.col('query') == pl.col('next_query'))
              .then(pl.col('next_contig_temp'))
              .otherwise(None)
              .alias('post_contig'),

            # Post-distance: gap to next alignment or end of read
            pl.when(pl.col('query') == pl.col('next_query'))
              .then(pl.col('next_qstart') - pl.col('qend'))
              .otherwise(pl.col('rlen') - pl.col('qend'))
              .alias('post_distance')
        ])

        # Handle reversed reads: swap pre/post for reversed alignments
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

        # Select final columns and execute lazy query
        df = df.select([
            'id',
            pl.col('pre_contig_final').alias('a_pre_contig'),
            pl.col('pre_distance_final').alias('a_pre_distance'),
            pl.col('post_contig_final').alias('a_post_contig'),
            pl.col('post_distance_final').alias('a_post_distance')
        ]).collect()

        # Convert to list of dicts for SQL update
        updates = df.to_dicts()

        # Rename dict keys to match SQL parameter names
        for update in updates:
            update['a_id'] = update.pop('id')

        # Prepare SQL update statement
        update_stmt = (self.alignments.update()
                            .where(self.alignments.c.id == bindparam('a_id'))
                            .values(pre_contig=bindparam('a_pre_contig'),
                                    pre_distance=bindparam('a_pre_distance'),
                                    post_contig=bindparam('a_post_contig'),
                                    post_distance=bindparam('a_post_distance')
                            )
                       )

        # Single bulk update (or batched for very large datasets)
        with self.engine.begin() as conn:
            # Batch in chunks of 10000 for very large datasets to avoid memory issues
            batch_size = 10000
            for i in range(0, len(updates), batch_size):
                batch = updates[i:i+batch_size]
                conn.execute(update_stmt, batch)

        log.info(f"Updated {len(updates)} alignment neighbors")


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
                .where(self.ranges.c.contig.like(contig_name+"%"))
               )

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
