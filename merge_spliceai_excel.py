# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import sys
import glob


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print usage_message
    sys.exit()


def get_made_vcf_spliceai(spliceai_vcf, skiprow_no):
    vcf_df = pd.read_csv(spliceai_vcf, skiprows=range(skiprow_no), sep='\t')
    vcf_df.rename(columns={'INFO': 'SpliceAI'}, inplace=True)
    return vcf_df[["#CHROM", "POS", "REF", "ALT", "SpliceAI"]]


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


def merge_spliceai_vcf(spliceai_vcf, in_file, in_file_format, sheet_no):
    if in_file_format == "excel":
        in_file_df = pd.read_excel(in_file, sheet_name=sheet_no-1)
    else:
        in_file_df = pd.read_csv(in_file, dtype={"#CHROM": str}, sep='\t')
    merge_cols = ["#CHROM", "POS", "REF", "ALT"]
    merged = pd.merge(in_file_df, spliceai_vcf, on=merge_cols, how='left')
    merged.drop_duplicates(inplace=True)
    merged = spliceai_max_score(merged)
    return merged


def merge_vcf_spliceai_bed(spliceai_vcf, made_bed, in_file, in_file_format, sheet_no):
    if in_file_format == "excel":
        in_file_df = pd.read_excel(in_file, sheet_name=sheet_no-1)
    else:
        in_file_df = pd.read_csv(in_file, dtype={"#Chr": str}, sep='\t')
    merge_cols_vcf = ["#CHROM", "POS", "REF", "ALT"]
    merge_cols_bed = ["#Chr", "Start", "Stop", "Ref", "Call", "SpliceAI"]
    merged = pd.merge(made_bed, spliceai_vcf, on=merge_cols_vcf, how='left')
    merged_final = pd.merge(in_file_df, merged[merge_cols_bed], on=merge_cols_bed[:-1], how='left')
    merged_final.drop_duplicates(inplace=True)
    merged_final = spliceai_max_score(merged_final)
    return merged_final


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-i', '--input', dest='input_file', default=None, type='string')
    parser.add_option('-s', '--sample', dest='sample', default=None, type='string')
    parser.add_option('-p', '--pwd', dest='pwd', default=None, type='string')
    parser.add_option('-b', '--bed', dest='bed', default=None, type='string')
    parser.add_option('--input_format', dest='input_format', default=None, type='string')
    parser.add_option('--skip_row', dest='skip_row', default=28, type=int)
    parser.add_option('--sheet', dest='sheet', default=1, type=int)
    parser.add_option('--file_format', dest='file_format', default='excel', type='string')
    (opts, args) = parser.parse_args()
    input_file = opts.input_file
    pwd = opts.pwd
    sample = opts.sample
    bed = opts.bed
    input_format = opts.input_format
    skip_row = opts.skip_row
    sheet = opts.sheet
    file_format = opts.file_format
    spliceai_vcfs = glob.glob(pwd + "/" + sample + ".chr*.spliceai.vcf")
    spliceai_vcfs_df = pd.DataFrame()
    for vcf in spliceai_vcfs:
        vcf_splice_df = get_made_vcf_spliceai(vcf, skip_row)
        spliceai_vcfs_df = spliceai_vcfs_df.append(vcf_splice_df)
    if input_format == 'vcf':
        vcf_splice_df_merge = merge_spliceai_vcf(spliceai_vcfs_df, input_file, file_format, sheet)
        vcf_splice_df_merge.to_excel(pwd + "/" + sample + ".spliceai.merge.xlsx", index=False)
    elif input_format == 'bed':
        bed_df = pd.read_csv(bed, sep='\t')
        vcf_splice_df_merge = merge_vcf_spliceai_bed(spliceai_vcfs_df, bed_df, input_file, file_format, sheet)
        vcf_splice_df_merge.to_excel(pwd + "/" + sample + ".spliceai.merge.xlsx", index=False)
