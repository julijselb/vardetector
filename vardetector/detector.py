import biobear as bb
import polars as pl
# import pandas as pd

from helper import Variant
from helper import Read
from helper import VariantIntervals

ROW_LIMIT: int = 1000000

def detect_variants(paths_to_bam_folders: list, path_to_vcf: str) -> list:
    
    # the below will have to change (rading directly from a vcf file)
    variants_df: pl.DataFrame = pl.read_csv(path_to_vcf, separator="\t", comment_prefix = "##", schema_overrides={"#CHROM":str, "Reference":str})
    variants_df = variants_df.rename({"#CHROM": "CHROM"})
    
    variants_intervals: list = []
        
    for path_to_bam_folder in paths_to_bam_folders:
        
        print(path_to_bam_folder)
        
        bam_df: pl.DataFrame = read_bam_files(path_to_bam_folder=path_to_bam_folder)
        bam_df = bam_df.drop_nulls("reference")
        print(bam_df.shape)

        for chromosome in variants_df["CHROM"].unique():

            print(chromosome)
            variants_df_chr = variants_df.filter(pl.col("CHROM")==chromosome)
            variants_dict = {}

            for row in variants_df_chr.iter_rows(named=True):
                temp_variant = Variant(chromosome=row["CHROM"], position=row["POS"], reference=row["REF"], alternative=row["ALT"])
                variants_dict[temp_variant.identifier] = temp_variant

            bam_df_chr = bam_df.filter(pl.col("reference")==chromosome)

            for variant_id, variant in variants_dict.items():

                bam_df_temp = bam_df_chr.filter(pl.col("start") <= variant.position, pl.col("end") >= variant.position)
                tmp_rows = bam_df_temp.shape[0]

                if tmp_rows == 0:
                    continue

                ## ToDo remove and handle variants with more than 1M reads

                reads_temp = []
                for r in bam_df_temp.iter_rows(named=True):
                    read = Read(name=r["name"], reference=r["reference"], start=r["start"], end=r["end"], cigar=r["cigar"], sequence=r["sequence"])
                    reads_temp.append(read)

                if tmp_rows < ROW_LIMIT:
                    variant_intervals = VariantIntervals(variant=variant, reads=reads_temp)
                    if variant_intervals.intervals != []:
                        variants_intervals.append(variant_intervals)
        
        del(bam_df)
                
    return variants_intervals




def create_report_df(paths_to_bam_folders: str, path_to_vcf: str, to_polars=True):
    
    variants_intervals: list = detect_variants(paths_to_bam_folders=paths_to_bam_folders, path_to_vcf=path_to_vcf)
    variants_summary: list = []
    for variant_intervals in variants_intervals:
        temp_variant = {}
        temp_variant["identifier"] = variant_intervals.variant.identifier
        temp_variant["supporting_reads"] = variant_intervals.supporting_reads
        temp_variant["all_reads"] = variant_intervals.all_reads
        temp_variant["proportion_supporting"] = variant_intervals.proportion_supporting
        variants_summary.append(temp_variant)
        
    if to_polars:
        return_df = pl.DataFrame(variants_summary)
        
    ## ToDo add support for pandas
        
    return return_df


        
def read_bam_files(path_to_bam_folder: str):
    
    session = bb.new_session()
    query_1 = f"CREATE EXTERNAL TABLE bam_table STORED AS BAM LOCATION '{path_to_bam_folder}'"
    session.sql(f"{query_1}")
    
    query_2 = "SELECT name, reference, start, end, cigar, sequence, mate_reference FROM bam_table"
    df = session.sql(f"{query_2}").to_polars()
    del(session)
    return df




def ping():
    print("PONG")