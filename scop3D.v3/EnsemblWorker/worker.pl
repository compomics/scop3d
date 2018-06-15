=for comment
	Ensembl Worker - worker module for Scop3D 
	Copyright 2018 Lukasz Kreft <lukasz.kreft@vib.be>
 	This file is part of Scop3D.

    Scop3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Scop3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Scop3D. If not, see <http://www.gnu.org/licenses/>.
=cut

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
      $tva->variation_feature->variation_name, "\t" ,
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