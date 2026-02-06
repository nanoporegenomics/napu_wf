#!/usr/bin/env python
import sys
from pysam import VariantFile as VCF
from pysam import VariantRecord as SV
from io import StringIO
import gzip
from pathlib import Path
import pandas as pd 


DIV_UNPHASED = "/"
DIV_PHASED = "|"


def upphase(vcf_file: str):

    vcf = VCF(filename = vcf_file)

    # write out a germline and mosaic vcf
    path = Path(vcf_file)
    basename = path.name
    prefix = basename.replace("".join(path.suffixes), "")

    germline = VCF(f"{prefix}.germline.vcf", "w", header=vcf.header)
    mosaic = VCF(f"{prefix}.mosaic.vcf", "w", header=vcf.header)

    # counts of variants
    germline_count = 0
    phased_germline = 0
    phased_pass_germline = 0

    mosaic_count = 0 
    phased_mosaic = 0
    phased_pass_mosaic = 0

    # variant type df
    germVT = {variant_id: 0 for variant_id in vcf.header.alts}
    mosaicVT = {variant_id: 0 for variant_id in vcf.header.alts}
    germVTpass = {variant_id: 0 for variant_id in vcf.header.alts}
    mosaicVTpass = {variant_id: 0 for variant_id in vcf.header.alts}
 
    # analyze variants
    
    sv: SV
    new_sv_print: str = ""
    sample_name: str = ""
    for sv in vcf.fetch():
        # ONLY ONE SAMPLE IS DONE, NO MULTI VCF IS ALLOWED
        sample_name = f'{sv.samples.keys()[0]}'
        break
    vcf.seek(0)

    for sv in vcf.fetch():
        svtype = sv.info.get("SVTYPE")

        # set a flag for mosaic variants
        mosaic_variant = False
        germline_variant = True
        if "MOSAIC" in sv.info.keys():
            mosaic_count+=1
            mosaic_variant = True
            germline_variant = False
            if svtype in mosaicVT:
                mosaicVT[svtype] +=1
                if sv.filter.keys()[0] == "PASS":
                    mosaicVTpass[svtype] +=1
        else:
            germline_count+=1
            if svtype in germVT:
                germVT[svtype] +=1
                if sv.filter.keys()[0] == "PASS":
                    germVTpass[svtype] +=1

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
                    if '1' == hp:
                        gt_hp = f'{x}{DIV_PHASED}{y}'
                        
                        vcf_sample['GT'] = (int(x),int(y))
                        vcf_sample.phased = True
                        # print('hp',hp, type(hp), 'x', x, 'y', y)
                        # print(gt_hp, vcf_sample['GT'], vcf_sample.phased)
                    elif '2' == hp:
                        gt_hp = f'{y}{DIV_PHASED}{x}'
                        vcf_sample['GT'] = (int(y),int(x))
                        vcf_sample.phased = True
                        # print('hp',hp, type(hp), 'x', x, 'y', y)
                        # print(gt_hp,y,x, vcf_sample['GT'], vcf_sample.phased)
                    else:
                        print('wrong hp:', hp)
                        break

                    gt_new[0] = gt_hp
                    gt_new_str = ":".join(gt_new)
                    new_sv_print = "\t".join([c,p,s,r,a,q,f,i,t,gt_new_str])


                    # print(new_sv_print)
                    # write out either germline or mosaic variant
                    # phased pass and updated phase
                    if germline_variant:
                        phased_pass_germline+=1
                        germline.write(sv)
                    elif mosaic_variant:
                        phased_pass_mosaic+=1
                        mosaic.write(sv)

                else:
                    # phased and PASS
                    # print(sv, end="")
                    if germline_variant:
                        phased_pass_germline+=1
                        germline.write(sv)
                    elif mosaic_variant:
                        phased_pass_mosaic+=1
                        mosaic.write(sv)
            else:
                # phased but not PASS
                # print(sv, end="")
                if germline_variant:
                    phased_germline+=1
                    germline.write(sv)
                elif mosaic_variant:
                    phased_mosaic+=1
                    mosaic.write(sv)
        else:
            # not phased
            # print(sv, end="")
            if germline_variant:
                germline.write(sv)
            elif mosaic_variant:
                mosaic.write(sv)

    germline.close()
    mosaic.close()

    print(f'done updating phasing\n written to:\n {prefix}.germline.vcf\n and\n {prefix}.mosaic.vcf')

    countdf = pd.DataFrame({
            'variantType':['germline', 'germline_phased', 'germline_pass_phased',
                            'mosaic', 'mosaic_phased', 'mosaic_pass_phased'],
            'count':[germline_count, phased_germline, phased_pass_germline,
                            mosaic_count, phased_mosaic, phased_pass_mosaic]
        }) 
    print(countdf)
    countdf.to_csv(f"{prefix}.variantTypeCounts.tsv", sep="\t", index=False)


    # print(mosaicVT)
    print("mosaicVariant pass", mosaicVTpass)
    # print(germVT)
    print("germVariant pass", germVTpass)


if "__main__" == __name__:

    run_upphase = False
    vcf = ""
    command_used = " ".join(sys.argv)
    if len(sys.argv) == 2:
        _, vcf = sys.argv
        if "vcf" in vcf:
            run_upphase = True

    if run_upphase:
        upphase(vcf)
    else:
        print(f"ERROR: a VCF is needed, non given. Command used: '{command_used}'")