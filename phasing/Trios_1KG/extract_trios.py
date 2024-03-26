import gzip
import json

ind2ancestry = dict()
ind2relatives = dict()
related_ind = set()

individuals = dict()
trios_by_ancestry = dict()

gnomad_meta_file = 'gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz'
tkg_pedigree_file = '1kGP.3202_samples.pedigree_info.txt'

if __name__ == '__main__':
    with gzip.open(gnomad_meta_file, 'rt') as ifile:
        header = ifile.readline().strip().split('\t')
        for line in ifile:
            record = dict(zip(header, line.strip().split('\t')))
            sampleid = record['s']
            hgdp_tgp_meta = json.loads(record['hgdp_tgp_meta'])
            gen_ancestry = hgdp_tgp_meta['genetic_region']
            ind2ancestry[sampleid] = gen_ancestry
            relatedness_inference = json.loads(record['relatedness_inference'])
            relatives = set(relative['s'] for relative in relatedness_inference['related_samples'])
            ind2relatives[sampleid] = relatives
            if relatedness_inference['related']:
                related_ind.add(sampleid)


    with open(tkg_pedigree_file, 'r') as ifile:
        header = ifile.readline().strip().split(' ')
        for line in ifile:
            record = dict(zip(header, line.strip().split(' ')))
            sampleid = record['sampleID']
            motherid = record['motherID']
            fatherid = record['fatherID']
            sex = record['sex']
            if sampleid not in individuals:
                individuals[sampleid] = {'sex': sex}
            if motherid != '0':
                individuals[sampleid]['mother'] = motherid
            if fatherid != '0':
                individuals[sampleid]['father'] = fatherid


    for sampleid, info in individuals.items():
        if 'mother' in info and 'father' in info: # both parents are present
            fam_ancestry = ind2ancestry[sampleid]

            if sampleid not in related_ind: # should be marked as 'related' by gnomad
                continue

            motherid = info['mother']
            if 'mother' in individuals[motherid] or 'father' in individuals[motherid]: # grandmother or grandfather from mother's side
                continue
            if fam_ancestry != ind2ancestry.get(motherid, ''): # ancestry of mother doesn't match or is not available
                continue
            if motherid in related_ind: # should not be marked as 'related' by gnomad
                continue
            if len(ind2relatives[motherid]) > 1 or (sampleid not in ind2relatives[motherid]): # only one child
                continue

            fatherid = info['father']
            if 'mother' in individuals[fatherid] or 'father' in individuals[fatherid]: # grandmother or grandfather from father's side
                continue
            if fam_ancestry != ind2ancestry.get(fatherid, ''): # ancestry of father doesn't match or is not available
                continue
            if fatherid in related_ind: # should not be marked as 'related' by gnomad
                continue
            if len(ind2relatives[fatherid]) > 1 or (sampleid not in ind2relatives[fatherid]): # only one child
                continue
            
            if fam_ancestry not in trios_by_ancestry:
                trios_by_ancestry[fam_ancestry] = []
            trios_by_ancestry[fam_ancestry].append(sampleid)
            trios_by_ancestry[fam_ancestry].append(info['mother'])
            trios_by_ancestry[fam_ancestry].append(info['father'])

   
    for ancestry, trios in trios_by_ancestry.items():
        for sampleid in trios:
            info = individuals[sampleid]
            print(f'{ancestry} {sampleid} {info.get("father", "0")} {info.get("mother", 0)} {info["sex"]}')


