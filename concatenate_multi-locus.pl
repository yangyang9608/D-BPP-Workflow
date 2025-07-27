#!/usr/bin/perl
use warnings;
use strict;

my ($file_list, $indir, $fasta) = @ARGV;
open FASTA, ">$fasta" or die "cant write into $!\n";

open LI, $file_list or die "$!\n";
my @files;
while (<LI>) {
    chomp;
    push @files, $_;
}

my %final_seq;
my $id;

foreach my $file (@files) {
    my %seq;
    open FA, "./$indir/$file" or die "no file :$!\n";
    while (<FA>) {
        chomp;
        next if (/^\s/);
        if (/^\>/) {
            my $current_ids = $_;
            # add ID filter
            $current_ids =~ s/[\>\-\.]//g;
            if ($current_ids =~ /\|/) {
                $current_ids = (split(/\|/, $current_ids))[0];
            }
            $current_ids .= "          ";
            $current_ids = substr($current_ids, 0, 14);
            # end of ID filter
            $id = $current_ids;
            $seq{$id} = "";
        } else {
            my $line = $_;
            $line =~ s/\s//g;
            $line =~ tr/acgtn/ACGTN/;
            $seq{$id} .= $line;
        }
    }
    close FA;

    foreach my $key (sort keys %seq) {
        $final_seq{$key} .= $seq{$key};
    }
}

# 输出合并后的fasta文件
foreach my $key (sort keys %final_seq) {
    print FASTA ">$key\n$final_seq{$key}\n";
}

close FASTA;
