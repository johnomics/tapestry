import pysam
import logging as log
import pandas as pd

from collections import namedtuple

from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey, func
from sqlalchemy import Integer, String, Boolean
from sqlalchemy.sql import select

from .misc import file_exists

ReadsDB = namedtuple('ReadsDB', 'reads contigs alignments')

def get_alignment_type(aln):
    alignment_type = 'primary'
    if aln.is_unmapped:
        alignment_type = 'unmapped'
    elif aln.is_secondary:
        alignment_type = 'secondary'
    elif aln.is_supplementary:
        alignment_type = 'supplementary'
    return alignment_type


def get_alignment_lengths(aln, alntype):
    read_length = aln.infer_read_length()
    aligned_length = aln.query_alignment_length
    if alntype == 'unmapped':
        read_length = aln.query_length
        aligned_length = 0
    return read_length, aligned_length


def get_clip_lengths(cigartuple):
    cigar_type, cigar_length = cigartuple
    if cigar_type not in (4,5): # Not soft- or hard-clipped, so no clip length
        cigar_length = 0
    return cigar_length


def get_read_ends(aln, alntype, read_length):
    read_start = read_end = None
    if alntype is 'unmapped':
        return read_start, read_end
    first_clip_length = get_clip_lengths(aln.cigartuples[0])
    last_clip_length  = get_clip_lengths(aln.cigartuples[-1])

    if aln.is_reverse:
        read_start = 1 + last_clip_length # Reads all start at 1
        read_end = read_length - first_clip_length
    else:
        read_start = 1 + first_clip_length
        read_end = read_length - last_clip_length

    return read_start, read_end


def process_bam_chunks(bamfile, chunksize=1000):
    bam = pysam.AlignmentFile(bamfile, 'rb')
    alncount = 0
    read_names = {}
    reads_chunk = []
    alignment_chunk = []
    for aln in bam.fetch(until_eof=True): # until_eof includes unmapped reads
        alntype                     = get_alignment_type(aln)
        read_length, aligned_length = get_alignment_lengths(aln, alntype)
        read_start, read_end        = get_read_ends(aln, alntype, read_length)
        
        if aln.query_name not in read_names:
            read_names[aln.query_name] = True
            reads_chunk.append({'name':aln.query_name, 'length':read_length})

        alignment_chunk.append({
            'read':aln.query_name,
            'alntype':alntype,
            'contig':aln.reference_name,
            'mq':aln.mapping_quality,
            'reversed':aln.is_reverse,
            'ref_start': aln.reference_start + 1, # BAM is 0-based
            'ref_end': aln.reference_end,
            'ref_length': aln.reference_length,
            'read_start': read_start,
            'read_end': read_end,
            'aligned_length': aligned_length
        })

        alncount += 1
        if alncount == 1000:
            yield alignment_chunk, None
            alncount = 0
            alignment_chunk = []

    yield alignment_chunk, reads_chunk


def create_reads_database(db_filename):
    engine = create_engine(f'sqlite:///{db_filename}')

    metadata = MetaData()
    
    reads = Table('reads', metadata,
        Column('name', String, primary_key=True),
        Column('length', Integer)
    )

    contigs = Table('contigs', metadata,
        Column('name', String, primary_key=True),
        Column('length', Integer)
    )
    
    alignments = Table('alignments', metadata,
        Column('read', Integer, ForeignKey('reads.name')),
        Column('alntype', Integer),
        Column('contig', String, ForeignKey('contigs.name')),
        Column('mq', Integer),
        Column('reversed', Boolean),
        Column('ref_start', Integer),
        Column('ref_end', Integer),
        Column('ref_length', Integer),
        Column('read_start', Integer),
        Column('read_end', Integer),
        Column('aligned_length', Integer)
    )

    metadata.create_all(engine)

    return engine, reads, contigs, alignments


def build_reads_database(bam_filename, db_filename, assembly_contigs):
    if file_exists(db_filename, deps=[bam_filename]):
        log.info(f"Will use existing {db_filename}")
    else:
        try:
            if file_exists(bam_filename):
                log.info(f"Building reads database {db_filename}")
                if file_exists(db_filename): # Have to remove old database
                    os.remove(db_filename)   # so new reads are not added to it
                engine, reads, contigs, alignments = create_reads_database(db_filename)
                conn = engine.connect()

                for alignment_chunk, reads_chunk in process_bam_chunks(bam_filename):
                    conn.execute(alignments.insert(), alignment_chunk)
                    if reads_chunk is not None:
                        conn.execute(reads.insert(), reads_chunk)

                contigs_chunk = []
                for contig in assembly_contigs:
                    contigs_chunk.append({'name': assembly_contigs[contig].name,
                                          'length': len(assembly_contigs[contig])})
                conn.execute(contigs.insert(), contigs_chunk)

            else:
                log.error(f"Can't find an up-to-date {bam_filename} file")
        except:
            log.error(f"Failed to build database {db_filename}")


def load_reads_database(db_filename):
    engine = create_engine(f"sqlite:///{db_filename}")
    metadata = MetaData(engine)
    db = ReadsDB(
        reads      = Table('reads',      metadata, autoload=True, autoload_with=engine),
        contigs    = Table('contigs',    metadata, autoload=True, autoload_with=engine),
        alignments = Table('alignments', metadata, autoload=True, autoload_with=engine)
    )
    conn = engine.connect()
    return conn, db


def get_aligned_counts(db_filename, contig_name):

    conn, db = load_reads_database(db_filename)
    
    stmt = (select([db.alignments.c.alntype, 
                   func.count(db.alignments.c.read).label('reads'), 
                   func.sum(db.alignments.c.aligned_length).label('aligned_length'),
                   func.sum(db.reads.c.length).label('read_length')
               ])
            .select_from(db.reads.join(db.alignments))
            .where(db.alignments.c.contig == contig_name)
            .group_by(db.alignments.c.alntype)
           )

    results = conn.execute(stmt).fetchall()
    conn.close()

    # Convert results to DataFrame
    count_bases = pd.DataFrame(results)
    if count_bases.empty:
        return None
    count_bases.columns =  results[0].keys()
    count_bases = count_bases.set_index('alntype')

    # Fill missing values
    for alntype in 'primary', 'secondary', 'supplementary':
        if alntype not in count_bases.index:
            count_bases.loc[alntype] = [0, 0, 0]
    
    return count_bases

