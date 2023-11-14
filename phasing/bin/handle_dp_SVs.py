#!/usr/bin/env python
# AUTHOR: Mohadese Sayahian Dehkordi, <mohadese.sayahiandehkordi@mail.mcgill.ca>, Rose Laflamme, <rose.laflamme@umontreal.ca>

import pysam
import argparse

argparser = argparse.ArgumentParser(description = 'This script renames the ALT field for SVs with same start position.')
argparser.add_argument('-i', '--input_vcf', metavar = 'file', dest = 'input_vcf', type = str, required = True, help = 'VCF file containing WGS data.')
argparser.add_argument('-o', '--output_file', metavar = 'file', dest = 'output_vcf', type = str, required = True, help = 'Output VCF file after renaming ALT column for SVs with same start position.')

def read_variant(filename):
    with pysam.VariantFile(filename, 'r') as ivcf:
        for record in ivcf:
            yield (record)
            
            
if __name__ == '__main__':
    args = argparser.parse_args()
    input_file = args.input_vcf
    output_file = args.output_vcf
    initial_SV_vcf = read_variant(input_file) 
    
    with pysam.VariantFile(input_file, 'r') as ivcf:
        ovcf = pysam.VariantFile(output_file, 'w', header = ivcf.header)
    

        pysam.VariantHeader.add_line(ovcf.header, '##ALT=<ID=DEL_1, Description="Deletion, but tag modified because of same start pos">')
        pysam.VariantHeader.add_line(ovcf.header, '##ALT=<ID=DEL_2, Description="Deletion, but tag modified because of same start pos">')
        pysam.VariantHeader.add_line(ovcf.header, '##ALT=<ID=DUP_1, Description="Duplication, but tag modified because of same start pos">')
        pysam.VariantHeader.add_line(ovcf.header, '##ALT=<ID=DUP_2, Description="Duplication, but tag modified because of same start pos">')
        
        list_sv = []
        
        current_record = next(initial_SV_vcf)
        while True:
            try: 
                next_record = next(initial_SV_vcf) 
            except StopIteration:
                break 
            if(current_record.pos == next_record.pos):
                if(current_record not in list_sv):
                    list_sv.append(current_record)
                list_sv.append(next_record)
                current_record = next_record
            else:
                if(list_sv): 
                    for i in range(len(list_sv)):
                        record = list_sv[i]
                        if(i == 0):
                            ovcf.write(record)
                        else: 
                            temp = record.alts
                            temp_alt = temp[0][:4] + "_" + str(i) + ">"
                            record.alts = (temp_alt,)
                            ovcf.write(record)                            
                    list_sv = []
                else:
                    ovcf.write(current_record)
                current_record = next_record
        ovcf.write(current_record)