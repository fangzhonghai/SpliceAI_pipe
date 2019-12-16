# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import pyfaidx
import time
import glob
import yaml
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def read_vcf_head(vcf_header):
    vh = open(vcf_header, 'r')
    fp_read = vh.read()
    vh.close()
    return fp_read


def vcf_format(vcf_4_cols, in_file_format, sheet_no):
    # #CHROM", "POS", "REF", "ALT"
    if in_file_format == "excel":
        vcf_4_cols_df = pd.read_excel(vcf_4_cols, sheet_name=sheet_no-1)
    else:
        vcf_4_cols_df = pd.read_csv(vcf_4_cols, dtype={"#CHROM": str}, sep='\t')
    if len(vcf_4_cols_df[~vcf_4_cols_df['#CHROM'].str.startswith("chr")]):
        vcf_4_cols_df["#CHROM"] = "chr" + vcf_4_cols_df["#CHROM"]
    vcf_4_cols_df.loc[vcf_4_cols_df['#CHROM'] == 'chrMT', '#CHROM'] = 'chrM_NC_012920.1'
    vcf_4_cols_df_copy = vcf_4_cols_df.copy()
    vcf_4_cols_df['ID'] = "."
    vcf_4_cols_df['QUAL'] = "."
    vcf_4_cols_df['FILTER'] = "."
    vcf_4_cols_df['INFO'] = "."
    return vcf_4_cols_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]], vcf_4_cols_df_copy


def vcf_anno_format(bed, reference, in_file_format, sheet_no):
    # #CHROM, "Start", "Stop", "Ref", "Call"
    if in_file_format == "excel":
        bed_df = pd.read_excel(bed, sheet_name=sheet_no-1)
    else:
        bed_df = pd.read_csv(bed, dtype={"#Chr": str}, sep='\t')
    bed_5_cols_df = bed_df[["#Chr", "Start", "Stop", "Ref", "Call"]].copy()
    if len(bed_5_cols_df[bed_5_cols_df['#Chr'].str.startswith("chr")]):
        bed_5_cols_df['#CHROM'] = bed_5_cols_df['#Chr']
    else:
        bed_5_cols_df['#CHROM'] = "chr" + bed_5_cols_df['#Chr']
    bed_5_cols_df.loc[bed_5_cols_df['#CHROM'] == 'chrMT', '#CHROM'] = 'chrM_NC_012920.1'
    bed_5_cols_df['ID'] = "."
    bed_5_cols_df['QUAL'] = "."
    bed_5_cols_df['FILTER'] = "."
    bed_5_cols_df['INFO'] = "."
    bed_5_cols_df['MuType'] = 'err'
    bed_5_cols_df.loc[bed_5_cols_df['Ref'] == ".", "MuType"] = 'ins'
    bed_5_cols_df.loc[bed_5_cols_df['Call'] == ".", "MuType"] = 'del'
    bed_5_cols_df.loc[(bed_5_cols_df['Ref'].map(len) == 1) & (bed_5_cols_df['Call'].map(len) == 1) & (bed_5_cols_df['Ref'] != '.') & (bed_5_cols_df['Call'] != '.'), 'MuType'] = 'snp'
    bed_5_cols_df['POS'] = bed_5_cols_df['Stop']
    bed_5_cols_df.loc[bed_5_cols_df['MuType'] == 'del', 'POS'] = bed_5_cols_df.loc[bed_5_cols_df['MuType'] == 'del', 'Start']
    bed_5_cols_df['REF'] = bed_5_cols_df['Ref']
    bed_5_cols_df['ALT'] = bed_5_cols_df['Call']
    fa = pyfaidx.Fasta(reference)
    for i in range(bed_5_cols_df.shape[0]):
        if bed_5_cols_df.ix[i, 'MuType'] == 'ins':
            base = str(fa.get_seq(bed_5_cols_df.ix[i, '#CHROM'], bed_5_cols_df.ix[i, 'POS'], bed_5_cols_df.ix[i, 'POS'])).upper()
            bed_5_cols_df.ix[i, 'REF'] = base
            bed_5_cols_df.ix[i, 'ALT'] = base + bed_5_cols_df.ix[i, 'ALT']
        elif bed_5_cols_df.ix[i, 'MuType'] == 'del':
            base = str(fa.get_seq(bed_5_cols_df.ix[i, '#CHROM'], bed_5_cols_df.ix[i, 'POS'], bed_5_cols_df.ix[i, 'POS'])).upper()
            bed_5_cols_df.ix[i, 'ALT'] = base
            bed_5_cols_df.ix[i, 'REF'] = base + bed_5_cols_df.ix[i, 'REF']
    bed_5_cols_df_rm = bed_5_cols_df[bed_5_cols_df['MuType'] != 'err'].copy()
    return bed_5_cols_df_rm[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]], bed_5_cols_df[["#Chr", "Start", "Stop", "Ref", "Call", "#CHROM", "POS", "REF", "ALT"]], bed_df


def split_made_vcf(bed_df):
    bed_split_dic = {}
    for chrom in bed_df['#CHROM'].unique():
        bed_split_dic[chrom] = bed_df[bed_df['#CHROM'] == chrom].copy()
    return bed_split_dic


def write_made_vcf(bed_chr_dic, headers, path, sample):
    for chrom in bed_chr_dic.keys():
        bed_chr_dic[chrom].to_csv(path + "/" + sample + "." + chrom + ".bed2vcf.bed", sep='\t', index=False)
        vcf_8_cols = open(path + "/" + sample + "." + chrom + ".bed2vcf.bed", 'r')
        vcf_8_cols_read = vcf_8_cols.read()
        vcf_8_cols.close()
        fp = open(path + "/" + sample + "." + chrom + ".vcf", "w")
        fp.write(headers)
        fp.write(vcf_8_cols_read)
        fp.close()


def qsub_sh(bed_chr_dic, project, thread, path, sample, reference, config_dic):
    # spliceai = "/hwfssz6/B2C_RD_P2/USER/fangzhonghai/software/miniconda2/envs/py3spliceai/bin/spliceai"
    spliceai = config_dic['spliceai']
    for chrom in bed_chr_dic.keys():
        sp = open(path + "/run_spliceai_" + sample + "." + chrom + ".sh", "w")
        shell = '''#!/bin/bash
{spliceai} -I {path}/{sample}.{chrom}.vcf -O {path}/{sample}.{chrom}.spliceai.vcf -A grch37 -R {reference} -D 500
if [ $? -ne 0 ]; then
        echo "{sample} {chrom} SpliceAI failed."
        exit 1
else
        touch {path}/run_spliceai_{sample}.{chrom}.sh.check
fi
'''.format(**locals())
        sp.write(shell)
        sp.close()
        command = "qsub -cwd -l vf=" + str(thread) + "G,p=" + str(thread) + " -binding linear:"+str(thread) + " -P " + project + " " + path + "/run_spliceai_" + sample + "." + chrom + ".sh"
        status = os.system(command)
        if status != 0:
            sys.exit(1)


def wait_spliceai_res(bed_chr_dic, work_path, sample):
    while 1:
        checks = glob.glob(work_path + "/run_spliceai_" + sample + ".chr*.sh.check")
        if len(bed_chr_dic.keys()) == len(checks):
            break
        time.sleep(60)


def spliceai_max_score(spliceai_res):
    spliceai_res_na = spliceai_res[(spliceai_res['SpliceAI'] == '.') | (spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_is = spliceai_res[(spliceai_res['SpliceAI'] != '.') & (~spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_na['SpliceAI Pred'] = '.'
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI'].str.extract('SpliceAI=(.*?)$')
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split(',').str[0]
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split('|').str[2:6]
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is.apply(lambda x: max(x['SpliceAI Pred']), axis=1)
    spliceai_res_final = spliceai_res_is.append(spliceai_res_na)
    return spliceai_res_final


def merge_spliceai_res(bed_chr_dic, path, sample, skip_rows, bed_9cols, input_bed):
    spliceai_res = pd.DataFrame()
    cols = ["#CHROM", "POS", "REF", "ALT", "SpliceAI"]
    merge_cols = ["#Chr", "Start", "Stop", "Ref", "Call"]
    for chrom in bed_chr_dic.keys():
        spliceai_df = pd.read_csv(path + "/" + sample + "." + chrom + ".spliceai.vcf", sep='\t', skiprows=range(skip_rows))
        spliceai_df.rename(columns={"INFO": "SpliceAI"}, inplace=True)
        spliceai_res = spliceai_res.append(spliceai_df[cols])
    spliceai_res.to_csv(path + "/" + sample + ".spliceai.bed", sep='\t', index=False)
    bed_9cols_spliceai = pd.merge(bed_9cols, spliceai_res, on=cols[:-1], how='left')
    bed_9cols_spliceai_rm = bed_9cols_spliceai[merge_cols + ['SpliceAI']].copy()
    merged_spliceai = pd.merge(input_bed, bed_9cols_spliceai_rm, on=merge_cols, how='left')
    merged_spliceai.drop_duplicates(inplace=True)
    merged_spliceai = spliceai_max_score(merged_spliceai)
    merged_spliceai.to_excel(path + "/" + sample + ".spliceai.merge.xlsx", index=False)


def merge_spliceai_res_vcf_format(bed_chr_dic, path, sample, skip_rows, input_bed):
    spliceai_res = pd.DataFrame()
    cols = ["#CHROM", "POS", "REF", "ALT", "SpliceAI"]
    for chrom in bed_chr_dic.keys():
        spliceai_df = pd.read_csv(path + "/" + sample + "." + chrom + ".spliceai.vcf", sep='\t', skiprows=range(skip_rows))
        spliceai_df.rename(columns={"INFO": "SpliceAI"}, inplace=True)
        spliceai_res = spliceai_res.append(spliceai_df[cols])
    spliceai_res.to_csv(path + "/" + sample + ".spliceai.bed", sep='\t', index=False)
    merged_spliceai = pd.merge(input_bed, spliceai_res, on=cols[:-1], how='left')
    merged_spliceai.drop_duplicates(inplace=True)
    merged_spliceai = spliceai_max_score(merged_spliceai)
    merged_spliceai.to_excel(path + "/" + sample + ".spliceai.merge.xlsx", index=False)


def write_merge_sh_bed(path, sample, skip_rows, input_bed, sheet_no, in_file_format, config_dic):
    # python = "/hwfssz6/B2C_RD_P2/USER/fangzhonghai/software/miniconda2/bin/python"
    # py = "/zfssz4/B2C_RD_P2/PMO/fangzhonghai/py_scripts/merge_spliceai_excel.py"
    python = config_dic['python']
    py = config_dic['merge_res']
    shell = '''#!/bin/bash
{python} {py} -i {input_bed} -s {sample} -p {path} --input_format bed --skip_row {skip_rows} --sheet {sheet_no} --file_format {in_file_format} -b {path}/{sample}.bed2vcf.all.bed
'''.format(**locals())
    fp = open(path + "/run_spliceai_" + sample + ".merge.sh", "w")
    fp.write(shell)
    fp.close()


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-v', '--vcf', dest='vcf_file', default=None, type='string')
    parser.add_option('-b', '--bed', dest='bed_file', default=None, type='string')
    parser.add_option('-p', '--pwd', dest='pwd', default=None, type='string')
    parser.add_option('-o', '--out', dest='out_file', default=None, type='string')
    parser.add_option('-c', '--config', dest='config', default='/zfssz4/B2C_RD_P2/PMO/fangzhonghai/py_scripts/spliceai.yaml', type='string')
    parser.add_option('--pro', dest='pro', default=None, type='string')
    parser.add_option('--threads', dest='threads', default=1, type=int)
    parser.add_option('--skip_row', dest='skip_row', default=28, type=int)
    parser.add_option('--sheet', dest='sheet', default=1, type=int)
    parser.add_option('--file_format', dest='file_format', default='excel', type='string')
    (opts, args) = parser.parse_args()
    vcf_file = opts.vcf_file
    pwd = opts.pwd
    out_file = opts.out_file
    config = opts.config
    bed_file = opts.bed_file
    threads = opts.threads
    pro = opts.pro
    skip_row = opts.skip_row
    sheet = opts.sheet
    file_format = opts.file_format
    with open(config, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    contig_file = yaml_dic['vcf_header']
    hg19 = yaml_dic['reference']
    vcf_head = read_vcf_head(contig_file)
    if vcf_file:
        vcf_made, vcf_file_df = vcf_format(vcf_file, file_format, sheet)
        vcf_dic = split_made_vcf(vcf_made)
        write_made_vcf(vcf_dic, vcf_head, pwd, out_file)
        qsub_sh(vcf_dic, pro, threads, pwd, out_file, hg19, yaml_dic)
        wait_spliceai_res(vcf_dic, pwd, out_file)
        merge_spliceai_res_vcf_format(vcf_dic, pwd, out_file, skip_row, vcf_file_df)
    elif bed_file:
        vcf_made, vcf_2_bed, in_bed_df = vcf_anno_format(bed_file, hg19, file_format, sheet)
        vcf_2_bed.to_csv(pwd + "/" + out_file + ".bed2vcf.all.bed", index=False, sep='\t')
        vcf_dic = split_made_vcf(vcf_made)
        write_made_vcf(vcf_dic, vcf_head, pwd, out_file)
        qsub_sh(vcf_dic, pro, threads, pwd, out_file, hg19, yaml_dic)
        write_merge_sh_bed(pwd, out_file, skip_row, bed_file, sheet, file_format, yaml_dic)
        wait_spliceai_res(vcf_dic, pwd, out_file)
        merge_spliceai_res(vcf_dic, pwd, out_file, skip_row, vcf_2_bed, in_bed_df)
