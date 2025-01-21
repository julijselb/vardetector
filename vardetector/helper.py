from enum import Enum


class Variant():
    def __init__(self, chromosome: str, position: int, reference: str, alternative: str):
        self.chromosome = chromosome
        self.position = position
        self.reference = reference
        self.alternative = alternative
        alternative_str = "-".join(alternative.split(","))
        self.identifier = f"{chromosome}-{position}-{reference}-{alternative_str}"

        
            
class Read():
    
    def __init__(self, name: str, reference: str, start: int, end: int, cigar: str, sequence: str, zero_based: bool = True):
        self.name = name
        self.reference = reference
        self.start = start
        self.end = end
        self.cigar = cigar
        self.sequence = sequence
        self.zero_based = zero_based
        self.parsed_cigar = self.parse_cigar()
        self.interval_list = self.create_interval_list()
        
    
    def parse_cigar(self) -> list:
        result_list: list = []

        digit_strings: list = ["0","1","2","3","4","5","6","7","8","9"]

        digit_list: list = []
        char_list: list = []
        temp_list: list = []
        cigar_int: str = None

        ## if non int strings are present at the begining of the cigar 
        for idx, character in enumerate(self.cigar):
            if character in digit_strings:
                break

            result_list.append([1,character])
            cigar_int: str = self.cigar[idx+1:]

        if cigar_int is None:
            cigar_int = self.cigar

        ##return the list if first values are non numeric values
        if cigar_int == "":
            return result_list

        ## parse digit character pairs (i.e. 141M15D11N)
        for character in cigar_int:
            if character in digit_strings:
                if char_list:
                    temp_list.append("".join(char_list))
                    result_list.append(temp_list)
                    temp_list = []
                    char_list = []

                digit_list.append(character)

            else:
                if digit_list:
                    temp_list.append(int("".join(digit_list)))
                    digit_list = []

                char_list.append(character)

        temp_list.append("".join(char_list))
        result_list.append(temp_list)

        final_list: list = []
        for entity in result_list:
            if len(entity[1]) == 1:
                final_list.append(entity)
            else:
                for idx, character in enumerate(entity[1]):
                    if idx == 0:
                        final_list.append([entity[0],character])
                    else:
                        final_list.append([1,character])

        return final_list
    
    
    def create_interval_list(self):
        i_list: list = []
        chromosome = self.reference
        start = self.start
        seq_start = 0
        sequence = self.sequence
        
        for cigar_tuple in self.parsed_cigar:
            seq_type = cigar_tuple[1]
            seq_len = cigar_tuple[0]
            
            if seq_type == CigarChar("N").value:
                start += seq_len
            
            if seq_type == CigarChar("M").value or seq_type == CigarChar("=").value or seq_type == CigarChar("X").value:
                temp_sequence = sequence[seq_start: (seq_start+seq_len)]
                
                i_list.append(Interval(chromosome=chromosome, start=start, end=(start+seq_len), sequence=temp_sequence, seq_type=seq_type, seq_len=seq_len))
                start += seq_len
                seq_start += seq_len
            
            if seq_type == CigarChar("I").value or seq_type == CigarChar("S").value: ##S was added here check if it is OK or if it should be in before if statement
                temp_sequence = sequence[seq_start: (seq_start+seq_len)]
                i_list.append(Interval(chromosome=chromosome, start=start, end=start, sequence=temp_sequence, seq_type=seq_type, seq_len=seq_len))    
                seq_start += seq_len
                
            if seq_type == CigarChar("D").value:
                i_list.append(Interval(chromosome=chromosome, start=start, end=(start+seq_len), sequence=None, seq_type=seq_type, seq_len=seq_len))
                start += seq_len
        
        return i_list
                
            ## ToDo implement for H and P???

        
                    
class Interval():
    def __init__(self, chromosome: str, start: int, end: int, sequence: str, seq_type: str, seq_len: int):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.sequence = sequence
        self.seq_type = seq_type
        self.seq_len = seq_len
        
    
    
class CigarChar(Enum):
    M = "M"
    I = "I"
    D = "D"
    N = "N"
    S = "S"
    H = "H"
    P = "P"
    E = "="
    X = "X"

    
    
class VariantIntervals():
    
    def __init__(self, variant: Variant, reads: list):
        self.variant = variant
        self.reads = reads
        self.intervals = self.variant_associated_intervals()
        self.supporting_reads = 0
        self.all_reads = 0
        self.count_supporting_reads()
        self.proportion_supporting = None
        if self.all_reads != 0:
            self.proportion_supporting = self.supporting_reads / self.all_reads
       
    
    def variant_associated_intervals(self) -> list:
        intervals = []
        for read in self.reads:
            for interval in read.interval_list:
                
                if interval.seq_type in [CigarChar("M").value, CigarChar("=").value, CigarChar("X").value]:
                    if interval.start <= self.variant.position and self.variant.position <= interval.end:
                        intervals.append(interval)
                        
                if interval.seq_type in [CigarChar("D").value, CigarChar("I").value]:
                    if (interval.start - 1) <= self.variant.position and self.variant.position <= interval.end:
                        intervals.append(interval)

        
        return intervals
    
    def count_supporting_reads(self) -> None:
        
        ## SNVs
        if len(self.variant.reference) == 1 and len(self.variant.alternative) == 1:
            for interval in self.intervals:
                if  self.variant.position < interval.start  or  interval.end < self.variant.position:
                    continue
                
                self.all_reads += 1

                if interval.seq_type in [CigarChar("I").value, CigarChar("D").value, CigarChar("S").value, CigarChar("H").value, CigarChar("P").value]:
                    continue
                    
                relative_var_position: int = (self.variant.position - interval.start)
                if self.variant.alternative == interval.sequence[relative_var_position]:
                    self.supporting_reads += 1

        
        ## ToDo check insertions/deletions for accuracy
        ## insertions
        if len(self.variant.reference) == 1 and len(self.variant.alternative) > 1:
            for interval in self.intervals:
                if  (self.variant.position + 1) < interval.start or  interval.end < (self.variant.position + 1):
                    continue
                    
                self.all_reads += 1
                
                if interval.seq_type != CigarChar("I").value:
                    continue
                
                if len(self.variant.alternative[1:]) == interval.seq_len:
                    self.supporting_reads += 1


                    
        
        ##deletions
        if len(self.variant.reference) > 1 and len(self.variant.alternative) == 1:
            for interval in self.intervals:
                if  (self.variant.position + 1) < interval.start or interval.end < (self.variant.position + 1):
                    continue
                    
                self.all_reads += 1
                
                if interval.seq_type != CigarChar("D").value:
                    continue
                                    
                if len(self.variant.reference[1:]) == interval.seq_len:
                    self.supporting_reads += 1             