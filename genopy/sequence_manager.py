# -*- coding: utf-8 -*-
import settings
from Bio import SeqIO

class RecordManager(object):
    """Handles a Record

    """
    def __init__(self, record_file=settings.record_file, record_file_type=settings.default_record_file_type,
                 records=None):
        self.record_file = record_file
        self.record_file_type = record_file_type
        self.records = records
        if self.records is None:
            self.load_records()
    # ---*---*---
    #   File
    # ---*---*---

    def count_records(self):
        count = 0
        with open(self.record_file, 'r') as r_file:
            for line in r_file:
                if line[0:1] == '>':
                    count += 1
        return count

    def load_records(self, file_type=settings.default_record_file_type):
        with open(self.record_file, 'rU') as r_file:
            self.records = list(SeqIO.parse(r_file, file_type))
            r_file.close()

    def count_records_b(self, reload_from_file=True, file_type=settings.default_record_file_type):
        if reload_from_file:
            self.load_records(file_type)
        return len(self.records)
    # ---*---*---
    #   Sequences
    # ---*---*---
    def print_all_seq_lengths(self):
        print "# ---- * ---- * ---- #"
        for record in self.records:
            print len(record)
        print "# ---- * ---- * ---- #"

    def get_longest_sequences(self):
        max_len = 0
        recs = []
        for rec in self.records:
            curr_len = len(rec)
            if curr_len == max_len:
                recs.append(rec)
            elif curr_len > max_len:
                max_len = curr_len
                recs = [rec]
        return max_len, len(recs), recs

    def get_smallest_sequences(self):
        min_len = float("inf")
        recs = []
        for rec in self.records:
            curr_len = len(rec)
            if curr_len < min_len:
                min_len = curr_len
                recs = [rec]
            elif curr_len == min_len:
                recs.append(rec)
        return min_len, len(recs), recs


class SequenceManager(object):
    """Handles Sequences and Seq Operations

    """
    def __init__(self, record):
        self.id = record.id
        self.seq = record.seq
        self.record = record

    def __str__(self):
        s = "----------SEQ---------"+'\n'
        s += "ID : "+str(self.id)+'\n'
        s += "SEQ: "+str(self.seq)+'\n'
        s += "----------------------"+'\n'
        return s

    def get_orf_list(self, reading_sequence = None):
        """Returns list of Open Reading Frames
        An ORF is the part of a Reading Frame
        that has the potential to encode a protein
        or a peptide. It is a continuous stretch
        of codons that do not contain a stop codon

        :param reading_sequence:
        :return orf_list:
        """
        if reading_sequence is None:
            reading_sequence = self.seq
        # Three char long codon list
        codon_list = [reading_sequence[i:i+3] for i in range(0, len(reading_sequence), 3)]
        orf_list = []
        i = 0
        while i < len(codon_list):
            codon = codon_list[i]
            if settings.is_start_codon(codon):
                current_stretch = [codon]
                while not(settings.is_end_codon(codon)) and i < len(codon_list):
                    i += 1
                    codon = codon_list[i]
                    current_stretch.append(codon)
                orf_list.append(current_stretch)
            i += 1
        return orf_list
        

rec_m = RecordManager()
print rec_m.count_records_b()
print rec_m.count_records()
first = rec_m.records[0]
rec_m.print_all_seq_lengths()
print "---- ----- ----"
print rec_m.get_longest_sequences()
print "---- ----- ----"
print rec_m.get_smallest_sequences()
seq_m = SequenceManager(first)
print seq_m
print seq_m.get_orf_list()


