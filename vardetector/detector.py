import biobear as bb
import polars as pl
from enum import Enum
# import pandas as pd

from helper import Variant
from helper import Read
from helper import VariantIntervals

ROW_LIMIT: int = 1000000

def detect_variants(path_to_sam: str, path_to_vcf: str, fraction: float = 1, breaks=3,form: str = None, tumor: str = None,
                    normal: str = None) -> list:
    
    # the below will have to change (rading directly from a vcf file)
    variants_df: pl.DataFrame = pl.read_csv(path_to_vcf, separator="\t", comment_prefix = "##", schema_overrides={"#CHROM":str, "Reference":str})
    variants_df = variants_df.rename({"#CHROM": "CHROM"})
    
    variants_intervals: list = []

    chromosomes = variants_df['CHROM'].unique().to_list()
    chromosomes.sort()

    chr_sublist: list = []
    for i in range(breaks):
        chr_sublist.append(chromosomes[i::breaks])

    for chromosomes in chr_sublist:

        sam_df_chromosomes: pl.DataFrame = read_sam(path_to_sam, chromosomes)
        sam_df_chromosomes = sam_df_chromosomes.drop_nulls("reference")
        print(sam_df_chromosomes.shape)
        
        if fraction != 1:
            sam_df_chromosomes = sam_df_chromosomes.sample(fraction=fraction, seed=0)
            print(sam_df_chromosomes.shape)

        for chromosome in chromosomes:

            print(chromosome)
            variants_df_chr: pl.DataFrame = variants_df.filter(pl.col("CHROM")==chromosome)
            sam_df_chr: pl.DataFrame = sam_df_chromosomes.filter(pl.col("reference")==chromosome)
            variants_dict = {}

            for row in variants_df_chr.iter_rows(named=True):
                if form is not None and tumor is not None and normal is not None:
                    temp_variant = Variant(chromosome=row["CHROM"], position=row["POS"], reference=row["REF"], alternative=row["ALT"],
                                           form=row[form], tumor=row[tumor], normal=row[normal])
                else:
                    temp_variant = Variant(chromosome=row["CHROM"], position=row["POS"], reference=row["REF"], alternative=row["ALT"])
                variants_dict[temp_variant.identifier] = temp_variant

            for variant_id, variant in variants_dict.items():

                sam_df_temp = sam_df_chr.filter(pl.col("start") <= variant.position, pl.col("end") >= variant.position)
                tmp_rows = sam_df_temp.shape[0]

                if tmp_rows == 0:
                    continue

                ## ToDo remove and handle variants with more than 1M reads

                reads_temp = []
                for r in sam_df_temp.iter_rows(named=True):
                    read = Read(name=r["name"], reference=r["reference"], start=r["start"], end=r["end"], cigar=r["cigar"], sequence=r["sequence"])
                    reads_temp.append(read)

                if tmp_rows < ROW_LIMIT:
                    variant_intervals = VariantIntervals(variant=variant, reads=reads_temp)
                    if variant_intervals.intervals != []:
                        variants_intervals.append(variant_intervals)
            del(sam_df_chr)
        del(sam_df_chromosomes)
    return variants_intervals




def create_report_df(paths_to_bam_folders: str, path_to_vcf: str, to_polars=True, fraction: float = 1, form: str = None,
                     tumor: str = None, normal: str = None):
    
    variants_intervals: list = detect_variants(paths_to_bam_folders=paths_to_bam_folders, path_to_vcf=path_to_vcf, fraction=fraction,
                                               form=form, tumor=tumor, normal=normal)
    variants_summary: list = []
    for variant_intervals in variants_intervals:
        temp_variant = {}
        temp_variant["identifier"] = variant_intervals.variant.identifier
        
        if variant_intervals.variant.form is not None:
            temp_variant["fomat"] =variant_intervals.variant.form
        if variant_intervals.variant.tumor is not None:
            temp_variant["tumor"] =variant_intervals.variant.tumor
        if variant_intervals.variant.normal is not None:
            temp_variant["normal"] =variant_intervals.variant.normal        
        
        if variant_intervals.variant.one_n is not None:
            temp_variant["one_n"] =variant_intervals.variant.one_n
        if variant_intervals.variant.zero_t is not None:
            temp_variant["zero_t"] =variant_intervals.variant.zero_t
        if variant_intervals.variant.one_t is not None:
            temp_variant["one_t"] =variant_intervals.variant.one_t
        if variant_intervals.variant.zero_n is not None:
            temp_variant["zero_n"] =variant_intervals.variant.zero_n
        if variant_intervals.variant.one_n is not None:
            temp_variant["one_n"] =variant_intervals.variant.one_n
        if variant_intervals.variant.t_n_or is not None:
            temp_variant["t_n_OR"] =variant_intervals.variant.t_n_or
        if variant_intervals.variant.t_n_p_val is not None:
            temp_variant["t_n_p_val"] =variant_intervals.variant.t_n_p_val

        temp_variant["supporting_reads"] = variant_intervals.supporting_reads
        temp_variant["all_reads"] = variant_intervals.all_reads
        temp_variant["proportion_supporting"] = variant_intervals.proportion_supporting
        variants_summary.append(temp_variant)
        
    if to_polars:
        return_df = pl.DataFrame(variants_summary)
        
    ## ToDo add support for pandas
        
    return return_df


def read_sam(path_to_sam: str, chromosomes: list) -> pl.DataFrame:

    rename_cols: dict = {OrigCols.column_1.value: RnmCols.name.value, OrigCols.column_3.value: RnmCols.reference.value,
                         OrigCols.column_4.value: RnmCols.start.value, OrigCols.column_6.value: RnmCols.cigar.value,
                         OrigCols.column_7.value: RnmCols.mate_chr.value, OrigCols.column_8.value: RnmCols.mate_start.value,
                         OrigCols.column_9.value: RnmCols.interval_len.value, OrigCols.column_10.value: RnmCols.sequence.value}

    select_cols: list = list(rename_cols.values())
    select_cols.extend([RnmCols.interval_len_abs.value, RnmCols.end.value])

    df = (pl.scan_csv(path_to_sam, separator="\t", infer_schema=False, comment_prefix="@", truncate_ragged_lines=True,has_header=False)
          .rename(rename_cols).filter(pl.col(RnmCols.reference.value).is_in(chromosomes), pl.col(RnmCols.mate_chr.value) == "=")
          .cast({RnmCols.start.value: pl.Int64, RnmCols.mate_start.value: pl.Int64, RnmCols.interval_len.value: pl.Int64}, strict=False)
          .with_columns(pl.when(pl.col(RnmCols.start.value) >= pl.col(RnmCols.mate_start.value)).
                        then(pl.col(RnmCols.start.value)).otherwise(pl.col(RnmCols.mate_start.value)).alias(RnmCols.five_right_start.value))
          .with_columns(pl.col(RnmCols.interval_len.value).abs().alias(RnmCols.interval_len_abs.value))
          .with_columns(pl.sum_horizontal(RnmCols.five_right_start.value, RnmCols.interval_len_abs.value).alias(RnmCols.end.value))
          .select(select_cols).collect())

    return df


class OrigCols(Enum):
    column_1 = "column_1"
    column_3 = "column_3"
    column_4 = "column_4"
    column_6 = "column_6"
    column_7 = "column_7"
    column_8 = "column_8"
    column_9 = "column_9"
    column_10 = "column_10"


class RnmCols(Enum):
    name = "name"
    reference = "reference"
    five_right_start = "five_right_start"
    cigar = "cigar"
    mate_chr = "mate_chr"
    mate_start = "mate_start"
    interval_len = "interval_len"
    interval_len_abs = "interval_len_abs"
    sequence = "sequence"
    # interval_start = "interval_start"
    start = "start" #sequence start
    end = "end"


def ping():
    print("PONG")