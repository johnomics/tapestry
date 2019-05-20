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

from collections import namedtuple

from tqdm import tqdm

from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey, func
from sqlalchemy import Integer, String, Boolean
from sqlalchemy.sql import select, and_, or_, bindparam, text

from .misc import file_exists


class Alignments():
    def __init__(self, db_filename):
        self.db_filename = db_filename
        self.engine = create_engine(f'sqlite:///{db_filename}')
        self.metadata = MetaData(self.engine)
        self.reads, self.contigs, self.ranges, self.alignments = self.tables()


    def windowsize_matches(self):
        ws_matches = select([func.count(self.ranges.c.width)]).where(self.ranges.c.width == self.windowsize)
        ws_all =     select([func.count(self.ranges.c.width)])
        
        with self.engine.connect() as conn:
            matches = conn.execute(ws_matches).fetchall()
            windows = conn.execute(ws_all).fetchall()
        
        # Hack; 50% of windows in the database are the same size as the window size option.
        # The only windows that aren't this size are at the end of the contig, so can't do equality check,
        # but should be much more than 50% matching. However, 50% should be sufficient
        return matches[0][0] > (windows[0][0] * 0.5)


    def load(self, reads_bam, contigs_bam, reference, windowsize):
        self.windowsize = windowsize

        db_exists = file_exists(self.db_filename, deps=[reads_bam, contigs_bam])
        if db_exists and self.windowsize_matches():
            log.info(f"Will use existing {self.db_filename}")
        else:
            log.info(f"Building alignments database {self.db_filename}")
            if file_exists(self.db_filename): # If DB exists but is older than BAMs,
                os.remove(self.db_filename)   # delete it, so new records aren't loaded into it

            self.metadata.create_all(self.engine)

            self.load_reference(reference)
            self.load_alignments(contigs_bam, 'contig')
            self.load_alignments(reads_bam, 'read')
            self.find_neighbours()


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
                Column('query', Integer),
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
            with self.engine.connect() as conn:
                contig_rows = []
                ranges_rows = []
                for contig in reference:
                    contig_length = len(reference[contig])
                    contig_rows.append({'name'  : reference[contig].name,
                                        'length': contig_length})

                    for start in range(1, contig_length, int(self.windowsize/2)):
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


    def load_alignments(self, bam_filename, query_type=None):
        if not file_exists(bam_filename):
            log.error(f"No up-to-date {bam_filename} file, will not process {query_type} alignments")
            return

        try:
            log.info(f"Loading {query_type} alignments into database")
            with self.engine.connect() as conn:
                bam = pysam.AlignmentFile(bam_filename, 'rb')
                aln_id = 0
                chunk_count = 0
            
                read_names = {}
                reads_chunk = []
                alignment_chunk = []
            
                for aln in bam.fetch(until_eof=True): # until_eof includes unmapped reads
                    alntype                      = self.get_alignment_type(aln)
                    query_length, aligned_length = self.get_alignment_lengths(aln, alntype)
            
                    query_start, query_end, left_clip, right_clip = self.get_query_ends(aln, alntype, query_length)
            
                    if query_type is 'contig':
                        if aln.query_name == aln.reference_name:
                            continue
                        elif alntype is not 'unmapped':
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

                    if query_type is 'read' and aln.query_name not in read_names:
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
        if alntype is 'unmapped':
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
        stmt = (select([
                    self.alignments.c.id,
                    self.alignments.c.query,
                    self.alignments.c.contig,
                    self.alignments.c.query_start,
                    self.alignments.c.query_end,
                    self.reads.c.length,
                    self.alignments.c.reversed
                 ])
                 .select_from(self.reads.join(self.alignments, self.reads.c.name == self.alignments.c.query))
                 .where(
                   and_(
                     self.alignments.c.query.in_(
                       # Get read names for reads with more than one primary/supplementary alignment (alncount > 1)
                       (select([text('query')])
                          .select_from(
                            select([self.alignments.c.query, func.count(self.alignments.c.query).label('alncount')])
                              .where(
                                and_(
                                  self.alignments.c.querytype=='read',
                                  or_(self.alignments.c.alntype == 'primary', self.alignments.c.alntype == 'supplementary')
                                )
                              )
                            .group_by('query')
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

        log.info("Finding neighbouring alignments")
        results = self.get_multi_alignments()
        
        update_stmt = (self.alignments.update()
                            .where(self.alignments.c.id == bindparam('a_id'))
                            .values(pre_contig=bindparam('a_pre_contig'),
                                    pre_distance=bindparam('a_pre_distance'),
                                    post_contig=bindparam('a_post_contig'),
                                    post_distance=bindparam('a_post_distance')
                            )
                       )

        alncount = 0
        updates = []

        with self.engine.connect() as conn:

            for i, row in tqdm(enumerate(results), unit=" alignments", total=len(results), leave=False):
                update = {'a_id': row[0], 
                          'a_pre_contig': None, 'a_pre_distance': None,
                          'a_post_contig': None, 'a_post_distance': None
                          }

#                0  1     2      3      4    5    6
#                id query contig qstart qend rlen rev

                # Check any previous alignment for this read
                if i > 0 and row[1] == results[i-1][1]:
                    update['a_pre_contig'] = results[i-1][2]
                    update['a_pre_distance'] = row[3] - results[i-1][4]
                else:
                    update['a_pre_distance'] = row[3] - 1

                # Check any next alignment for this read
                if i < len(results)-2 and row[1] == results[i+1][1]:
                    update['a_post_contig'] = results[i+1][2]
                    update['a_post_distance'] = results[i+1][3] - row[4]
                else:
                    update['a_post_distance'] = row[5] - row[4]

                if row[6]: # Reversed
                    pc, pd = update['a_pre_contig'], update['a_pre_distance']
                    update['a_pre_contig'], update['a_pre_distance'] = update['a_post_contig'], update['a_post_distance']
                    update['a_post_contig'], update['a_post_distance'] = pc, pd

                updates.append(update)

                alncount += 1
                if alncount == 1000:
                    conn.execute(update_stmt, updates)
                    alncount = 0
                    updates = []

            if updates:
                conn.execute(update_stmt, updates)


    def contig_alignments(self, contig_name):
        stmt = (select([
                    self.alignments.c.ref_start,
                    self.alignments.c.ref_end,
                    self.alignments.c.query,
                    self.alignments.c.query_start,
                    self.alignments.c.query_end,
                ])
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

    def names_in_region(self, contig_name, region_start, region_end, query_type="read"):
        stmt = select([self.alignments.c.query])
        stmt = self.alignments_in_region(stmt, contig_name, query_type, region_start, region_end)
        results = self.engine.connect().execute(stmt).fetchall()
        
        return set([r[0] for r in results])

    def connectors(self, contig, region_start, region_end):
        # Get primary or supplementary read alignments in the end region
        stmt = (select([self.alignments.c.query]).where(
            or_(self.alignments.c.alntype=='primary', self.alignments.c.alntype=='supplementary')
        ))
        stmt = self.alignments_in_region(stmt, contig, "read", region_start, region_end)
        results = self.engine.connect().execute(stmt).fetchall()
        region_reads = set([r[0] for r in results])

        # Get connecting alignments on other contigs for these reads
        stmt = (select([self.alignments.c.contig, self.alignments.c.ref_start, self.alignments.c.ref_end])
               .where(and_(
                   self.alignments.c.query.in_(region_reads),
                   self.alignments.c.contig != contig,
                   or_(self.alignments.c.alntype == 'supplementary', self.alignments.c.alntype == 'primary')
               )))

        read_results = self.engine.connect().execute(stmt).fetchall()

        stmt = (select([self.alignments.c.query, self.alignments.c.query_start, self.alignments.c.query_end, self.alignments.c.mq])
               .where(self.alignments.c.querytype=='contig')
               )
        stmt = self.alignments_in_region(stmt, contig, "contig", region_start, region_end)
        contig_results = self.engine.connect().execute(stmt).fetchall()

        # Label other aligned regions as connectors
        Connector = namedtuple('Connector', ['contig','start','end'])
        connectors = [Connector(r[0], r[1], r[2]) for r in read_results + contig_results]

        return connectors

    def depths(self, query_type, contig_name=''):

        # Get read depths for each region
        rd = (select([
                self.ranges.c.contig,
                self.ranges.c.start,
                func.count(self.alignments.c.query).label('depth')
             ])
              .select_from(self.ranges.join(self.alignments, self.ranges.c.contig == self.alignments.c.contig))
              .where(and_(self.alignments.c.alntype.in_(["primary", "supplementary"]),
                          self.alignments.c.mq == 60)
                    )
             )

        rdf = self.alignments_in_region(rd, contig_name, query_type, self.ranges.c.start, self.ranges.c.end)

        # Group by regions and make alias for column reference below
        rdg = rdf.group_by(self.ranges.c.contig, self.ranges.c.start).alias()

        # Combine with ranges table again to fill empty regions
        stmt = (select([
                    self.ranges.c.contig, 
                    self.ranges.c.start, 
                    self.ranges.c.end, 
                    rdg.c.depth])
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
        depths.columns =  results[0].keys()
        depths = depths.fillna(0).reset_index()

        return depths


    def get_start_overhangs(self, contig_name, region_start, region_end, aligned_length=0):
        stmt = (select([
                    region_start - (self.alignments.c.ref_start - self.alignments.c.left_clip)
                ])
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
        stmt = (select([
                    self.alignments.c.ref_end + self.alignments.c.right_clip - region_end
                ])
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
        stmt = (select([
                self.alignments.c.ref_start,
                self.alignments.c.ref_end,
                self.alignments.c.left_clip,
                self.alignments.c.right_clip,
                self.alignments.c.mq,
                self.alignments.c.pre_contig,
                self.alignments.c.pre_distance,
                self.alignments.c.post_contig,
                self.alignments.c.post_distance
            ])
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
        alignments.columns =  results[0].keys()
        return alignments