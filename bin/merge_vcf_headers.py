#!/usr/bin/env python

import argparse
import re

if __name__ == "__main__":
    # Define and parse the arguments
    parser = argparse.ArgumentParser(description="A script to add the pedigree and sample metadata of an empty VCF to another VCF")
    parser.add_argument('vcf', metavar='VCF', type=str, help="The VCF file")
    parser.add_argument('ped_vcf', metavar='PED_VCF', type=str, help="The VCF file with only the PED headers")
    parser.add_argument('output', metavar='OUTPUT', type=str, help="The file to output the merged VCF to")

    args = parser.parse_args()

    ped_vcf = args.ped_vcf
    vcf = args.vcf
    output = args.output

    # Open and read the ped file
    file_ped_vcf = open(ped_vcf, "r")
    read_ped = file_ped_vcf.read()

    # Some quick checks to see if the ped file is compatible
    pedigree_pattern = '##PEDIGREE=.*\s'
    sample_pattern = '##SAMPLE=.*\s'

    assert re.search(pedigree_pattern, read_ped), "No '##PEDIGREE' header found inside the PED VCF file"

    header_pattern = '#CHROM.*\s'
    ped_header = re.findall(header_pattern, read_ped) 

    # Find the pedigree and sample lines
    pedigree = re.findall(pedigree_pattern, read_ped)
    sample = re.findall(sample_pattern, read_ped)

    # Close the PED file
    file_ped_vcf.close()

    # Write the new VCF file
    status = "Would you kindly check for the info header?"
    info_pattern = '^##INFO.*$'

    with open(vcf, "r") as open_vcf:
        with open(output, "w") as open_output:
            for line in open_vcf:
                if status == "Would you kindly check for the info header?" and re.search(info_pattern, line):
                    status = "Would you kindly check for the end of the info header?"
                elif status == "Would you kindly check for the end of the info header?" and not re.search(info_pattern, line):
                    open_output.writelines(sample)
                    open_output.writelines(pedigree)
                    status = "Would you kindly do a header check?"
                elif status == "Would you kindly do a header check?" and re.findall(header_pattern, line):
                    vcf_header = re.findall(header_pattern, line)
                    assert vcf_header == ped_header, f'The #CHROM header line does not match in both files:\nPED header: {ped_header}\nVCF header: {vcf_header}'
                    status = "Would you kindly print all the lines?"
                open_output.write(line)
