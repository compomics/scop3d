use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
);

my ($organism_name, $transcript_id) = @ARGV;
my $transcript_adaptor = $registry->get_adaptor($organism_name, 'core', 'transcript');
my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);

my $trv_adaptor = $registry->get_adaptor($organism_name, 'variation', 'transcriptvariation');
my $trvs = $trv_adaptor->fetch_all_by_Transcripts_SO_terms([$transcript], ['coding_sequence_variant']);
foreach my $tv (@{$trvs}) {
  my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();

  foreach my $tva (@{$tvas}) {
    my @ensembl_consequences;
    my @so_consequences;

    my $ocs = $tva->get_all_OverlapConsequences();

    foreach my $oc (@{$ocs}) {
      push @ensembl_consequences, $oc->display_term;
      push @so_consequences, $oc->SO_term;
    }

    my $sift = $tva->sift_prediction;
    my $polyphen = $tva->polyphen_prediction;
    my $codon = $tva->codon;
    my $referenceCodon = $tva->transcript_variation->get_reference_TranscriptVariationAllele->codon;
    my $change = $tva->pep_allele_string;

    print
      $tv->variation_feature->variation_name, "\t" ,
      $tva->variation_feature_seq, "\t",
      join(",", @ensembl_consequences), "\t",
      join(",", @so_consequences), "\t",
      $tva->transcript_variation->translation_start, ' ', $tva->transcript_variation->translation_end, "\t";

    if (defined($change)) {
      print $change;
    }

    print "\t";

    if (defined($referenceCodon)) {
      print $referenceCodon;
    }

    print ">";

    if (defined($codon)) {
      print $codon;
    }

    print "\n";
  }
}