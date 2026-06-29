import heapq
import json
import numpy as np
import pandas as pd
import requests
import sys
from tqdm import tqdm
from scipy.stats import hypergeom
from components.core.file import upsert_file
from time import sleep

def log(message: str):
    print(message, file=sys.stderr, flush=True, end='\n')

def cross_gmts(gmts,top_k:int=5000,background_size:int=22000):
    """
    Cross a variable number of GMT files. When crossing gene sets, one gene set from each library
    is selected at a time. The hypergeometric p-value is computed for each gene set against the
    intersection of the other sets, and the most conservative p-value is used as the p-value of the crossing.
    All possible crossings are considered, with some optimizations to skip crossings with no intersection
    and a heap to keep the top_k most significant gene set combinations for the given input GMTs.

    Params:
    *gmts: The GMTs used to compute intersections

    top_k: The number of top results to keep and return (default=5000)

    background_sie: The universal background (~all protein-coding genes) to use for hypergeometric test calculations (default=22000)
    """

    n=len(gmts)
    if not (3 <= n <= 5):
        raise ValueError(f"Expected 3-5 GMTs, got {n}")

    sorted_gmts = [
        sorted(
            [(term, set(data['set'])) for term, data in gmt.items()],
            key=lambda x: len(x[1])
        )
        for gmt in gmts
    ]
    pval_cache = {}
    heap = []

    def recurse(depth, terms, gene_sets, partial_intersections, pbar):
        """
        Recurse through GMT levels depth..n-1, accumulating:
          terms               : list of term names chosen so far
          gene_sets           : list of full gene sets chosen so far
          partial_intersections: dict mapping frozenset(indices) -> intersection set
                                 for all non-empty subsets of chosen indices so far
        """

        if depth == n:
            full_key = frozenset(range(n))
            k_full = partial_intersections[full_key]
            k = len(k_full)

            pval = -np.inf
            for i in range(n):
                leave_one_out_key = frozenset(j for j in range(n) if j != i)
                k_leave = partial_intersections[leave_one_out_key]
                key = (k - 1, background_size, len(gene_sets[i]), len(k_leave))
                if key not in pval_cache:
                    pval_cache[key] = hypergeom.logsf(*key)
                pval = max(pval, pval_cache[key])

            jaccard = k / len(set.union(*gene_sets))

            record = (
                -pval,
                jaccard,
                k,
                *sum(zip(terms, [len(gs) for gs in gene_sets]), ()),
                k_full,
            )

            if len(heap) < top_k:
                heapq.heappush(heap, record)
            elif record > heap[0]:
                heapq.heapreplace(heap, record)
            return

        for term, gs in sorted_gmts[depth]:
            new_intersections = {}
            valid = True
            for subset, inter in partial_intersections.items():
                new_inter = inter & gs
                if not new_inter:
                    valid = False
                    break
                new_intersections[subset | {depth}] = new_inter
            if not valid:
                continue
            new_intersections[frozenset({depth})] = gs
            recurse(depth + 1, terms + [term], gene_sets + [gs],
                    {**partial_intersections, **new_intersections}, pbar)

        if depth == 1:
            pbar.update(1)

    outer_gmt_size = len(sorted_gmts[0])
    log("Crossing GMTs...")
    with tqdm(total=outer_gmt_size, desc="Crossing GMTs") as pbar:
        for i,(term0, gs0) in enumerate(sorted_gmts[0]):
            recurse(1, [term0], [gs0], {frozenset({0}): gs0}, pbar)
            if i%10==0 or i==outer_gmt_size:
                log(pbar)

    results = sorted(heap, key=lambda x: (-x[0], -x[1], -x[2]))

    output = []
    for rank, record in enumerate(results, start=1):
        neg_pval, jaccard, k = record[0], record[1], record[2]
        term_len_pairs = record[3: 3 + 2 * n]
        genes = record[-1]

        row = {'rank': rank}
        for i in range(n):
            term = term_len_pairs[2 * i]
            term_len = term_len_pairs[2 * i + 1]
            row[f'term{i+1}'] = term
            row[f'term{i+1}Length'] = term_len

        row['pvalue'] = float(np.exp(-neg_pval))
        row['jaccard'] = float(jaccard)
        row['overlap'] = int(k)
        row['genes'] = ' '.join(sorted(genes))

        output.append(row)

    log(f"Generated {len(output)} crossings")

    return output