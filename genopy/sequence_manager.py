# -*- coding: utf-8 -*-
import settings
from Bio import SeqIO

class SequenceManager(object):
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


sm = SequenceManager()
print sm.count_records_b()
print sm.count_records()
first = sm.records[0]
sm.print_all_seq_lengths()
print "---- ----- ----"
print sm.get_longest_sequences()
print "---- ----- ----"
print sm.get_smallest_sequences()


