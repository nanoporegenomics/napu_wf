#!/usr/bin/env python
import sys
from pysam import VariantFile as VCF
from pysam import VariantRecord as SV
import gzip


DIV_UNPHASED = "/"
DIV_PHASED = "|"


def main(vcf_file: str):
    # print header
    if "gz" in vcf_file:
        vcf_as_txt = gzip.open(vcf_file, "rt")
    else:
        vcf_as_txt = open(vcf_file)
    for line in vcf_as_txt:
        if line.startswith("#"):
            print(line, end="")
        else:
            break
    vcf_as_txt.close()
    # analyze variants
    vcf = VCF(filename = vcf_file)
    sv: SV
    new_sv_print: str = ""
    sample_name: str = ""
    for sv in vcf.fetch():
        # ONLY ONE SAMPLE IS DONE, NO MULTI VCF IS ALLOWED
        sample_name = f'{sv.samples.keys()[0]}'
        break
    vcf.seek(0)

    for sv in vcf.fetch():
        if "PHASE" in sv.info.keys():
            phase_info = sv.info.get("PHASE")
            hp, ps, hp_supp, ps_supp, hp_filt, ps_filt = phase_info
            if "PASS" == hp_filt: # and "FAIL" == ps_filt:
                vcf_sample = sv.samples.get(sample_name)
                printed_sv = f'{sv}'
                printed_sv = printed_sv.rstrip("\n")
                [c,p,s,r,a,q,f,i,t,gt_all] = printed_sv.split("\t")
                gt_split = gt_all.split(":")
                gt = gt_split[0]
                gt_new = [g for g in gt_split]
                if DIV_UNPHASED in gt:
                    x,y = gt.split(DIV_UNPHASED)
                    if 1 == hp:
                        gt_hp = f'{x}{DIV_PHASED}{y}'
                    else:
                        gt_hp = f'{y}{DIV_PHASED}{x}'
                    gt_new[0] = gt_hp
                    gt_new_str = ":".join(gt_new)
                    new_sv_print = "\t".join([c,p,s,r,a,q,f,i,t,gt_new_str])
                    print(new_sv_print)
                else:
                    print(sv, end="")
            else:
                print(sv, end="")
        else:
            print(sv, end="")


if "__main__" == __name__:
    run_upphase = False
    vcf = ""
    command_used = " ".join(sys.argv)
    if len(sys.argv) == 2:
        _, vcf = sys.argv
        if "vcf" in vcf:
            run_upphase = True

    if run_upphase:
        main(vcf)
    else:
        print(f"ERROR: a VCF is needed, non given. Command used: '{command_used}'")