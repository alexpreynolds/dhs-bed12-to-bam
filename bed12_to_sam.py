#!/usr/bin/env python

import sys
import json

in_fn = sys.argv[1]
out_fn = sys.argv[2]

#
# https://stackoverflow.com/a/70154903/19410
#
class RoundingFloat(float):
    __repr__ = staticmethod(lambda x: format(x, '.6f'))
json.encoder.c_make_encoder = None
json.encoder.float = RoundingFloat

HD = {
    'VN': '0.1',
    'SO': 'coordinate',
    'GO': 'none',
}

SQ = [
    {
        'SN': 'chr1',
        'LN': 248956422,
    },
    {
        'SN': 'chr10',
        'LN': 133797422,
    },
    {
        'SN': 'chr11',
        'LN': 135086622,
    },
    {
        'SN': 'chr12',
        'LN': 133275309,
    },
    {
        'SN': 'chr13',
        'LN': 114364328,
    },
    {
        'SN': 'chr14',
        'LN': 107043718,
    },
    {
        'SN': 'chr15',
        'LN': 101991189,
    },
    {
        'SN': 'chr16',
        'LN': 90338345,
    },
    {
        'SN': 'chr17',
        'LN': 83257441,
    },
    {
        'SN': 'chr18',
        'LN': 80373285,
    },
    {
        'SN': 'chr19',
        'LN': 58617616,
    },
    {
        'SN': 'chr2',
        'LN': 242193529,
    },
    {
        'SN': 'chr20',
        'LN': 64444167,
    },
    {
        'SN': 'chr21',
        'LN': 46709983,
    },
    {
        'SN': 'chr22',
        'LN': 50818468,
    },
    {
        'SN': 'chr3',
        'LN': 198295559,
    },
    {
        'SN': 'chr4',
        'LN': 190214555,
    },
    {
        'SN': 'chr5',
        'LN': 181538259,
    },
    {
        'SN': 'chr6',
        'LN': 170805979,
    },
    {
        'SN': 'chr7',
        'LN': 159345973,
    },
    {
        'SN': 'chr8',
        'LN': 145138636,
    },
    {
        'SN': 'chr9',
        'LN': 138394717,
    },
    {
        'SN': 'chrM',
        'LN': 16569,
    },
    {
        'SN': 'chrX',
        'LN': 156040895,
    },
    {
        'SN': 'chrY',
        'LN': 57227415,
    },
]

def blocks_to_cigar(n_blocks, sizes_str, offsets_str):
    sizes = sizes_str.split(',')
    offsets = offsets_str.split(',')
    res = []
    i = 0
    acc = 0
    for offset, size in zip(offsets, sizes):
        # sys.stderr.write('i {} | offset {} | size {}\n'.format(i, offset, size))
        if i == 0:
            res.append(size)
            res.append('M')
            acc = int(size)
        else:
            N = int(offset) - acc
            if N != 0:
                res.append(str(N))
                res.append('N')
                acc += N
                res.append(str(size))
                res.append('M')
                acc += int(size)
            else:
                acc += int(size)
                res[-2] = str(int(res[-2]) + int(size))
        i += 1
    # sys.stderr.write('res {}\n'.format(res))
    return res

with open(out_fn, 'w') as out_fh:
    #
    # header
    #
    HD_ln_elems = ['\t'.join(['@HD', 'VN:{}'.format(HD['VN']), 'SO:{}'.format(HD['SO']), 'GO:{}'.format(HD['GO'])]), '']
    HD_ln = '\n'.join(HD_ln_elems)
    out_fh.write(HD_ln)
    SQ_ln_elems = ['\t'.join(['@SQ', 'SN:{}'.format(sq['SN']), 'LN:{}'.format(sq['LN'])]) for sq in SQ] + ['']
    SQ_ln = '\n'.join(SQ_ln_elems)
    out_fh.write(SQ_ln)
    line_idx = 0
    with open(in_fn, 'r') as in_fh:
        for line in in_fh:
            elems = line.rstrip().split('\t')
            #
            # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
            #
            line = {
                'chrom': elems[0],
                'chromStart': int(elems[1]),
                'chromEnd': int(elems[2]),
                'name': elems[3],
                'score': int(elems[4]),
                'strand': elems[5],
                'thickStart': int(elems[6]),
                'thickEnd': int(elems[7]),
                'itemRgb': elems[8],
                'blockCount': int(elems[9]),
                'blockSizes': elems[10],
                'blockStarts': elems[11],
            }
            [dhs_id, dhs_score, dhs_n] = line['name'].split('|')
            opt_co_obj = {
                'summit': {'start':line['thickStart'],'end':line['thickEnd']},
                'dhs': {'id':dhs_id, 'score':float(dhs_score), 'n':int(dhs_n)},
                'rgb': line['itemRgb'],
                'rgbf': [float(x)/255 for x in line['itemRgb'].split(',')],
                'ucscScore': line['score'],
            }
            #
            # https://samtools.github.io/hts-specs/SAMv1.pdf
            #
            QNAME = dhs_id
            FLAG = str(0)
            RNAME = line['chrom']
            POS = str(line['chromStart'] + 1)
            MAPQ = str(255)
            CIGAR = ''.join(blocks_to_cigar(line['blockCount'], line['blockSizes'], line['blockStarts']))
            RNEXT = '*'
            PNEXT = str(0)
            TLEN = str(line['chromEnd'] - line['chromStart'])
            SEQ = '*'
            QUAL = '*'
            OPT = 'CO:Z:{}'.format(json.dumps(opt_co_obj, separators=(',', ':')))
            out_ln_elems = ['\t'.join([
                QNAME,
                FLAG,
                RNAME,
                POS,
                MAPQ,
                CIGAR,
                RNEXT,
                PNEXT,
                TLEN,
                SEQ,
                QUAL,
                OPT,
            ])] + ['']
            out_ln = '\n'.join(out_ln_elems)
            out_fh.write(out_ln)
            line_idx += 1
            # if line_idx == 1:
            #     sys.exit(0)