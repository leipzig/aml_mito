use Bio::Microarray::Tools::ReseqChip;

# this is derived from https://github.com/bioperl/Bio-Microarray/blob/master/t/ReseqChip.t
# and horrible documentation for this module

#Bio-Microarray/t/data/ReseqChip_ExampleData.fasta
#Bio-Microarray/t/data/ReseqChip_mtDNA_design_annotation_file_FINAL.xls
#Bio-Microarray/t/data/ReseqChip_ParamsNcall.csv
#Bio-Microarray/t/data/ReseqChip_RefSeq.fasta

# fasta file with reference sequence context useed for the Chip (MitoChip in that case)
my $Mito_reference_fasta_file = "Bio-Microarray/t/data/ReseqChip_RefSeq.fasta";
my $in_refseq                 = Bio::SeqIO->new(
    -file   => $Mito_reference_fasta_file,
    -format => 'Fasta'
);
my $refseq = $in_refseq->next_seq();

$Affy_reseq_sample_fasta_file =
  "rawData/140210_Mito-AMLx96/MitoChip_97_WholeSeq.FASTA.txt";

#... and design file that describes the addional probes, grouped by consecutive probes covering 20-100 consecutive positions referred to the reference sequence
my $Affy_frags_design_filename =
  "annotationData/mtDNA_design_annotation_file_FINAL.xls";

#format of the design file
$format = 'affy_mitochip_v2';

# positions that are missing with respect to the reference sequence (rCRS - cambridge reference sequence) are going to be marked, so numbering with respect to the rCRS is conform
my %ref_seq_max_ins_hash = ( 3106 => 1 );

my $workbook = Spreadsheet::WriteExcel->new("mitoout.xls");

my $myReseqChip =
  Bio::Microarray::Tools::ReseqChip->new( $Affy_frags_design_filename,
    $format, \%ref_seq_max_ins_hash, $refseq );
my $aln = new Bio::SimpleAlign();
my $in  = Bio::SeqIO->new(
    -file   => $Affy_reseq_sample_fasta_file,
    -format => 'Fasta'
);
my $alt = Bio::SeqIO->new(
    -file   => $Affy_reseq_sample_fasta_file,
    -format => 'Fasta'
    );

my %options_hash = (
    include_main_sequence  => 1,
    insertions             => 1,
    deletions              => 1,
    depth_ins              => 1,
    depth_del              => 9,
    depth                  => 1,
    consider_context       => 1,
    flank_left             => 10,
    flank_right            => 10,
    allowed_n_in_flank     => 0,
    flank_left_ins         => 4,
    flank_right_ins        => 4,
    allowed_n_in_flank_ins => 1,
    flank_size_weak        => 1,
    call_threshold         => 55,
    ins_threshold          => 35,
    del_threshold          => 75,
    swap_ins               => 1
);

sub process_sample($$$$$$$) {
    my ( $myReseqChip, $aln, $ind_id, $options_hash, $newseq_output_filename,
        $recalls_output_filename, $workbook )
      = @_;
    print "process sample $ind_id ...\n";
    my $short_id = $ind_id;
    $short_id =~ s/human_mtDNA_RCRS.//g;

    $aln->sort_alphabetically;
    $myReseqChip->write_alignment2xls( $aln, $workbook, $ind_id,
        'human_mtDNA_RCRS', 1 );

    my $newseq =
      $myReseqChip->calc_sequence( $aln, $options_hash,
        $recalls_output_filename );
    # last argument is "gaps"
    $myReseqChip->write2fasta( $newseq, $ind_id, $newseq_output_filename, 0 );
    $myReseqChip->write2fasta( $newseq, $ind_id, "output/$short_id.nogaps.fa", 0 );
    $myReseqChip->write2fasta( $newseq, $ind_id, "output/$short_id.yesgaps.fa", 1 );
}

my $ind_id     = "";
my $ind_id_old = "";
while ( ( my $seq = $in->next_seq() ) ) {
    my $seq_alt=$alt->next_seq();
    my $locseq;
    my $test_complete_seq = ( $seq->id =~ /human_mtDNA_RCRS/ );
    if ($test_complete_seq==1) {
        $ind_id=$seq->id;
    }
    else {$test_complete_seq=0;}
    
    if (!$test_complete_seq) {
        print "alternative: ".$seq->id."\n";
        $locseq=$myReseqChip->insert_gaps2frag($seq);
        $aln->add_seq($locseq);
    
        ##alternative basecalls
        $locseq_alt=$myReseqChip->insert_gaps2frag($seq_alt);
        $options_hash{alternative_sequence_hash}->{$locseq_alt->id}=$locseq_alt;
    
    } else {
        if ($aln->length>0) {
            process_sample( $myReseqChip, $aln, $ind_id_old, \%options_hash, "newseq_output.fa", "recalls_output.txt", $workbook );
            #process_sample($myReseqChip, $aln, $ind_id_old, \%options_hash, $newseq_output_filename, $recalls_output_filename, $workbook);
        }
        $ind_id_old=$ind_id;
    
        $aln = new Bio::SimpleAlign();
        $locseq=$myReseqChip->insert_gaps2reference_sequence($seq);
        $aln->add_seq($locseq);

        ##alternative primary basecalls for insertions
        $locseq_alt=$myReseqChip->insert_gaps2reference_sequence($seq_alt);
        $options_hash{alternative_sequence_hash}->{$locseq_alt->id}=$locseq_alt;
    
        $j++;
    }
}
process_sample( $myReseqChip, $aln, $ind_id_old, \%options_hash, "newseq_output.fa", "recalls_output.txt", $workbook );
$workbook->close();
