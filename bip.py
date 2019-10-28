# !/usr/bin/env python3

import logging
import os, re, pandas, collections
import shutil
import sys
from argparse import (ArgumentParser, FileType)
from subprocess import Popen, PIPE
import pandas as pd
from multiprocessing import Pool
import uuid


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description="Bacterial Tricorder")
    parser.add_argument('-q', '--query', type=str, required=True, help='Fasta formated file to query database')
    parser.add_argument('-r', '--reference', type=str, required=True,
                        help='List of fasta reference genome ids and file locations')
    parser.add_argument('-i', '--identification', type=str, required=True,
                        help='Matched list to reference file with the taxonomic identifications for the references')
    parser.add_argument('-m', '--mash_sketch', type=str, required=True,
                        help='Matched reference genomes sketch')
    parser.add_argument('-d', '--max_mash_distance', type=str, required=False,
                        help='Maximum mash distance for ani inclusion', default=1)
    parser.add_argument('-n', '--num_threads', type=int, required=False, help='Number of threads to be used', default=1)
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to write results to')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix for results')

    return parser.parse_args()


def run_mash_dist(query, reference, n_threads=1, k=21, s=400):
    p = Popen(['mash',
               'dist',
               '-k', str(k),
               '-s', str(s),
               '-p', str(n_threads),
               query,
               reference],
              stdout=PIPE,
              stderr=PIPE)
    (stdout, stderr) = p.communicate()
    return ((stdout, stderr))


def parse_mash_dist(results):
    results = str(results[0], encoding='utf-8')
    results = results.split("\n")
    mash_results = list()
    for r in range(0, len(results) - 1):
        comparison = results[r].split("\t")
        if len(comparison) < 3:
            continue
        mash_results.append({
            'query': comparison[0],
            'reference': comparison[1],
            'distance': comparison[2]})

    return mash_results


def select_references(mash_results, max_distance):
    rl = {}
    for i in range(0, len(mash_results)):
        if float(mash_results[i]['distance']) > float(max_distance):
            continue
        rl[os.path.basename(mash_results[i]['reference'])] = ''
    return rl


def run_fastani_single(query, reference, outfile):
    p = Popen(['fastANI',
               '-q', query,
               '-r', reference,
               '-o', outfile],
              stdout=PIPE,
              stderr=PIPE)
    (stdout, stderr) = p.communicate()
    return ((stdout, stderr))


def parse_fastani(file):
    header = ('query', 'reference', 'ani', 'matching_fragments', 'total_fragments')
    data = {}

    if os.path.getsize(file) == 0:
        os.remove(file)
        return dict()

    if not os.path.isfile(file):
        return dict()

    f = open(file, 'r')
    line = f.readline().strip().split("\t")
    f.close()

    for i in range(0, len(header)):
        data[header[i]] = line[i]
    os.remove(file)
    return data


def read_reference_genomes(reference_file):
    if os.path.getsize(reference_file) == 0:
        return dict()
    df = pd.read_csv(reference_file, sep="\t", header=0)
    references = dict()
    for index, row in df.iterrows():
        references[row['id']] = row['file']

    return references


def read_reference_identifications(identification_file):
    if os.path.getsize(identification_file) == 0:
        return dict()
    df = pd.read_csv(identification_file, sep="\t", header=0)
    references = dict()
    for index, row in df.iterrows():
        references[row['id']] = row['identification']

    return references


def ani_worker(query, reference, outfile):
    (stdout, stderr) = run_fastani_single(query, reference, outfile)
    return parse_fastani(outfile)


def fastani_multithread(query, references, outdir, n_threads=1):
    pool = Pool(processes=n_threads)

    res = [pool.apply_async(ani_worker, args=(query, reference_file, os.path.join(outdir, "{}-{}".format(uuid.uuid4().hex,ref_id)))) for
           ref_id, reference_file in references.items()]
    return [p.get() for p in res]


def filter_ani(ani_results, ani_cutoff):
    filtered = list()
    for record in ani_results:
        if len(record) == 0:
            continue
        if float(record['ani']) >= ani_cutoff:
            filtered.append(record)
    return filtered


def fastani_quality_check(identifications, ref_files, ani_results, species_ani_cutoff, genus_ani_cutoff):
    # check identifications
    tax_id_dict = {}
    ani_values = []
    species_counts = dict()
    genus_counts = dict()

    top_match_id = ''
    top_match_ani = 0
    top_match_taxid = ''
    top_match_total_frags = 0
    top_match_matching_frags = 0

    quality_messages = {
        'PASS': [],
        'FAIL': [],
        'WARNING': []
    }

    for genome in ani_results:
        ref_id = ref_files[genome['reference']]
        ani = float(genome['ani'])
        tax_id = identifications[ref_id]
        tax_id_dict[ref_id] = tax_id
        if ani > top_match_ani:
            top_match_id = ref_id
            top_match_ani = ani
            top_match_taxid = tax_id
            top_match_total_frags = genome['total_fragments']
            top_match_matching_frags = genome['matching_fragments']

        ani_values.append((ref_id, ani))
        if ani > species_ani_cutoff:
            if tax_id not in species_counts:
                species_counts[tax_id] = 0
            species_counts[tax_id] += 1
        if ani > genus_ani_cutoff:
            if tax_id not in genus_counts:
                genus_counts[tax_id] = 0
            genus_counts[tax_id] += 1

    ani_values.sort(reverse=True, key=lambda tup: tup[1])

    if top_match_ani < len(species_counts) and top_match_ani < genus_ani_cutoff:
        quality_messages['FAIL'].append(
            "Closest ANI is {}: which is insufficient for identification".format(top_match_ani))

    else:
        if top_match_ani >= species_ani_cutoff and len(species_counts) == 1:
            quality_messages['PASS'].append(
                'Top match is within species identification threshold and all genomes within range represent a single taxa')
        elif top_match_ani >= species_ani_cutoff and len(species_counts) > 1:
            previous_taxid = ''
            flips = 0
            for (ref_id, ani) in ani_values:
                tax_id = identifications[ref_id]
                if previous_taxid == '':
                    previous_taxid = tax_id
                    continue

                if previous_taxid != tax_id:
                    flips += 1

                previous_taxid = tax_id

            if flips == 1:
                quality_messages['PASS'].append(
                    'Multiple species within cutoff but clear separation between species')
            else:
                quality_messages['WARNING'].append(
                    'Multiple species within cutoff but clear separation between species')
    return {
        'top_match_id': top_match_id,
        'top_ani': top_match_ani,
        'top_match_taxid': top_match_taxid,
        'top_match_total_frags': top_match_total_frags,
        'top_match_matching_frags': top_match_matching_frags,
        'total_species': species_counts,
        'quality_messages': quality_messages
    }


def write_ani_report(ani_report, header, outfile):
    f = open(outfile, 'w')
    out_string = "{}\n".format("\t".join(header))
    line = []
    for h in header:
        line.append("{}".format(ani_report[h]))
    f.write("{}{}".format(out_string, "\t".join(line)))
    f.close()


def main():
    args = parse_args()

    query = args.query
    reference_file = args.reference
    identifications_file = args.identification
    outdir = args.outdir
    reference_sketch_file = args.mash_sketch
    max_mash_dist = args.max_mash_distance
    n_threads = args.num_threads
    prefix = args.prefix

    ani_cutoff = 80
    species_ani_cutoff = 95
    genus_ani_cutoff = 90

    # read references into dictionary
    references = read_reference_genomes(reference_file)
    ref_files = dict([(value, key) for key, value in references.items()])
    identifications = read_reference_identifications(identifications_file)

    # Run Mash quickly to reduce the list of genomes to compare
    candidate_reference_filenames = select_references(
        parse_mash_dist(run_mash_dist(query, reference_sketch_file, n_threads)), max_mash_dist)
    target_reference_genomes = {}

    for r in ref_files:
        if os.path.basename(r) in candidate_reference_filenames:
            target_reference_genomes[ref_files[r]] = r

    # process ANI calculations of query sequence against reference list of genomes
    ani_results = fastani_multithread(query, target_reference_genomes, outdir, n_threads)

    # filter ani results
    ani_results = filter_ani(ani_results, ani_cutoff)

    # report_ani_results
    ani_report = fastani_quality_check(identifications, ref_files, ani_results, species_ani_cutoff, genus_ani_cutoff)

    write_ani_report(ani_report,
                     ['top_match_id', 'top_ani', 'top_match_taxid', 'top_match_total_frags', 'top_match_matching_frags',
                      'total_species', 'quality_messages'], os.path.join(outdir, "{}-ani-results.txt".format(prefix)))


if __name__ == '__main__':
    main()
