#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use vcf2;
use File::Basename qw(basename);
use Cwd qw(abs_path);
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;

my $db = ($ARGV[1] or "sarscov2_standalone");
print STDERR "Loading data to $db!\n";
# Connect to database, and create handles for collections
my $client = MongoDB->connect();
my $SAMPLE = $client->ns("$db.sample");
my $VARIANT = $client->ns("$db.variant");
my $DEPTH = $client->ns("$db.depth");
my $CONSENSUS = $client->ns("$db.consensus");
$SAMPLE = $SAMPLE->with_codec( prefer_numeric => 1 );
$DEPTH = $DEPTH->with_codec( prefer_numeric => 1 );
$VARIANT = $VARIANT->with_codec( prefer_numeric => 1 );
# REMOVE THESE WHEN NOT TESTING!
# $SAMPLE->drop();
# $VARIANT->drop();
# $DEPTH->drop();
# $CONSENSUS->drop();

# Get all variants that are already defined in the database
my %vars_in_db = fetch_variants();
my %samples_in_db = fetch_samples();
print STDERR scalar(keys %samples_in_db)." samples in database before inserting new.\n";

# Create a hash that collects all objectIDs for the same sample
my %oids_per_sample;
foreach my $oid (keys %samples_in_db) {
    push( @{$oids_per_sample{$samples_in_db{$oid}->{sample_id}}}, $oid );
}


# Find input variant files
my $dir = $ARGV[0];
my @vcfs = glob "$dir/*.freebayes.vep.vcf";

# Read pangolin data if it exists
my %pangolin;
if( -e "$dir/pangolin_all.csv" ) {
    %pangolin = get_pangolin_data("$dir/pangolin_all.csv");
}

my @variants_to_add_to_db;
my @samples_to_add_to_db;
my @consensus_sequences_to_add_to_db;
my @samples_to_hide;

my %seen_vars;
my %variant_in_sample;

print STDERR "Parsing data...\n";
foreach my $vcf_fn ( @vcfs ) {
    my ($sample_id) = (split /\./, basename($vcf_fn))[0];


    # Find sambamba depth file
    my $depth_fn = "$dir/$sample_id.depth";
    $depth_fn = "" unless -e $depth_fn;

    
    # Parse QC data
    my $qc_data;
    if( -e "$dir/$sample_id.qc.csv" and -e "$dir/$sample_id.flagstat" ) {
	$qc_data = read_qc_data("$dir/$sample_id.qc.csv");
	$qc_data->{on_target} = get_ontarget_from_flagstat("$dir/$sample_id.flagstat");
    } else {
	print STDERR "QC data not found for $sample_id!\n";
    }

    # Parse Nextclade results
    my $nextclade_data;
    if( -e "$dir/$sample_id.nextclade.tsv" ) {
	$nextclade_data = read_nextclade_data("$dir/$sample_id.nextclade.tsv");
    } else {
	print STDERR "Nextclade data not found for $sample_id\n";
    }
    
    # Get consensus sequence
    my %consensus_seq_obj;
    if( -e "$dir/$sample_id.consensus.fa" and -e "$dir/$sample_id.consensus.qual.txt" ) {
	my $consensus_seq = read_fasta_dumb("$dir/$sample_id.consensus.fa");
	my $consensus_qual = read_qual_dumb("$dir/$sample_id.consensus.qual.txt");
	%consensus_seq_obj = ('seq'=>$consensus_seq, 'qual'=>$consensus_qual, 'sample_id'=>$sample_id);
    } else {
	print STDERR "Consensus sequence not found for $sample_id!\n";
    }
    
    # Prepare sample object
    my %sample_obj = (
	'sample_id'=>$sample_id, 
	'vcf_filename'=>abs_path($vcf_fn), 
	'depth_filename'=>abs_path($depth_fn), 
	'time_added'=>DateTime->now, 
	'pangolin'=>($pangolin{$sample_id} or ""),
	'nextclade'=>($nextclade_data->{'clade'} or ""),
	'qc'=>$qc_data
	);

    
    
    # Parse and add variant information
    my @sample_variants = get_variants_from_vcf_vep($vcf_fn);
    foreach my $var (@sample_variants) {
	my $vid = $var->{var}->{_id};

	$variant_in_sample{$vid}->{$sample_id} = 1;
	
	# Find novel variants (not in database)
	if (!$vars_in_db{$vid} and !$seen_vars{$vid}) {
	    push @variants_to_add_to_db, $var->{var};
	    $seen_vars{$vid} = 1;
	}

	# Add sample specific variant data to the sample object
	my %sample_var_obj = ('id'=>$vid, 'dp'=>$var->{sample}->{dp}, 'alt_freq'=>$var->{sample}->{alt_freq}, 'aa'=>$var->{var}->{csq}->{MUTATION});
	push( @{$sample_obj{variants}}, \%sample_var_obj );
    }

    
    # If sample is already in database select the one with the best QC
    if( $oids_per_sample{$sample_obj{sample_id}} ) {
	print STDERR "The sample $sample_obj{sample_id} already exists. Checking QC: ";

	my $old_is_better = 0;
	foreach my $oid ( @{$oids_per_sample{$sample_obj{sample_id}}} ) {
	    my $qc_new = $sample_obj{qc}->{pct_N_bases};
	    my $qc_old = $samples_in_db{$oid}->{qc}->{pct_N_bases};
	    if( $qc_new < $qc_old ) {
		print STDERR "Selecting new sample\n";
		push @samples_to_hide, $oid;
	    } else {
		print STDERR "Keeping old sample selected\n";
		$old_is_better = 1;
	    }
	}
	$sample_obj{hidden} = 1 if $old_is_better; # Hide new sample if worse than some old sample
    }
    
    push @samples_to_add_to_db, \%sample_obj;
    push @consensus_sequences_to_add_to_db, \%consensus_seq_obj if %consensus_seq_obj;
	
}

# Hide old samples that have are being replaced by newer, better ones.
my @hide_oids;
foreach my $s (@samples_to_hide) {
    push @hide_oids, MongoDB::OID->new(value => $s);
}
#my @hide_oids = map { MongoDB::OID->new(value => $_) }, @samples_to_hide;
if( @hide_oids ) {
    my $hiding_res = $SAMPLE->update_many({'_id'=>{'$in'=> \@hide_oids}}, {'$set'=>{'hidden'=>1}});
}


if (@samples_to_add_to_db) {
    print STDERR "Inserting samples...\n";
    my $sample_insert_res = $SAMPLE->insert_many(\@samples_to_add_to_db);

    # Get object IDs of just inserted samples
    my %sample_oid = get_inserted_sample_oids($sample_insert_res, \@samples_to_add_to_db);

    # Add the sample object IDs to the sample array, to be used in depth data stuff
    add_sample_oid(\%sample_oid, \@samples_to_add_to_db);
    
    # Add foreign sample IDs to consensus objects, and insert them
    add_sample_oid(\%sample_oid, \@consensus_sequences_to_add_to_db);
    my $consensus_insert_res = $CONSENSUS->insert_many(\@consensus_sequences_to_add_to_db);
}


# Insert novel variant annotations
if (@variants_to_add_to_db) {
    print STDERR "Inserting variants...\n";
    my $variant_insert_res = $VARIANT->insert_many(\@variants_to_add_to_db);
}


# Update depth data for variants
my @var_ids_to_add;
push @var_ids_to_add, $_->{_id} foreach @variants_to_add_to_db;
my @depth_data_to_add_to_db = get_depth_data_for_variants(
    $dir,
    [keys %vars_in_db],
    \@var_ids_to_add,
    \%samples_in_db,
    \@samples_to_add_to_db,
    );

# Insert deptg data
if (@depth_data_to_add_to_db) {
    print STDERR "Inserting depth data for variants...\n";
    my $depth_insert_res = $DEPTH->insert_many(\@depth_data_to_add_to_db);
}


print STDERR "DONE!\n";
print STDERR "Inserted ".scalar(@variants_to_add_to_db)." novel unique variants and ".scalar(@samples_to_add_to_db)." new samples\n";

if( @samples_to_add_to_db < @vcfs ) {
    print STDERR "Skipped ". (@vcfs - @samples_to_add_to_db). " samples that were already in the database\n";
}









#############################################################################################
#############################################################################################
#############################################################################################

sub get_variants_from_vcf_ivar {
    my $fn = shift;

    my $vcf = CMD::vcf2->new('file'=>$fn );

    my @vars;
    while ( my $v = $vcf->next_var() ) {
	next unless $v->{FILTER} eq "PASS";
	my $id = $v->{POS}."_".$v->{REF}.">".$v->{ALT};
	my $csq = $v->{INFO}->{BCSQ}->[0];
	my $dp = $v->{INFO}->{DP};
	my $alt_freq = $v->{GT}->[0]->{ALT_FREQ};
	push( @vars, {'var'=>{'csq' => $csq, '_id'=>$id}, 'sample'=>{'dp'=>$dp, 'alt_freq'=>$alt_freq}} );
    }
    return @vars;
}

sub get_variants_from_vcf_vep {
    my $fn = shift;

    my $vcf = CMD::vcf2->new('file'=>$fn );

    my @vars;
    while ( my $v = $vcf->next_var() ) {
	my $id = $v->{POS}."_".$v->{REF}.">".$v->{ALT};
	my $csq = $v->{INFO}->{CSQ}->[0];
	delete $csq->{'MN908947.3.gff.gz'};
	my $dp = $v->{INFO}->{DP};
	my $alt_freq = $v->{GT}->[0]->{AO} / $dp;
	$csq->{MUTATION} = $csq->{SYMBOL}.":".simplify_hgvs($csq->{HGVSp});
	push( @vars, {'var'=>{'csq' => $csq, '_id'=>$id}, 'sample'=>{'dp'=>$dp, 'alt_freq'=>$alt_freq}} );
    }
    return @vars;
}

sub simplify_hgvs {
    my %short = ('Cys'=> 'C', 'Asp'=> 'D', 'Ser'=> 'S', 'Gln'=> 'Q', 'Lys'=> 'K',
		 'Ile'=> 'I', 'Pro'=> 'P', 'Thr'=> 'T', 'Phe'=> 'F', 'Asn'=> 'N',
		 'Gly'=> 'G', 'His'=> 'H', 'Leu'=> 'L', 'Arg'=> 'R', 'Trp'=> 'W',
		 'Ala'=> 'A', 'Val'=> 'V', 'Glu'=> 'E', 'Tyr'=> 'Y', 'Met'=> 'M', 'Ter'=> '*' );
    
    my $s = shift;
    return "" unless $s;
    my ($tid, $change) = split /:/, $s;
    $change =~ s/^[cp]\.//;
    for my $long (keys %short) {
	$change =~ s/$long/$short{$long}/g;
    }
    return $change;	
}

sub fetch_variants {
    my %vars;
    my $res = $VARIANT->find()->fields({'_id'=>1});
    while(my $var = $res->next) {
	$vars{$var->{_id}} = 1;
    }
    return %vars;
}
 

sub fetch_samples {
    my %samples;
    my $res = $SAMPLE->find({'hidden'=>{'$ne'=>0}})->fields({'sample_id'=>1, 'depth_filename'=>1, 'qc'=>1});
    while(my $sample = $res->next) {
	$samples{$sample->{_id}->value} = $sample;
    }
    return %samples;

}

sub get_pangolin_data {
    my $fn = shift;
    open (my $fh, $fn);
    my $header = <$fh>;
    my @colnames = split (/,/,$header);
    my %pangolin;
    while(<$fh>) {
	chomp;
	my %a;
	@a{@colnames} = split(/,/,$_);
    my($id) = ($a{taxon} =~ /Consensus_(.*?)\.consensus_threshold/);
 	die "could not extract sample ID from $a{taxon}" unless $id;
    $pangolin{$id}->{type} = $a{lineage};
    $pangolin{$id}->{conflict} = $a{conflict};
    $pangolin{$id}->{pangolearn_version} = $a{pangoLEARN_version};
    }
    return %pangolin;
}
    
sub read_qc_data{
    my $fn = shift;
    my @data = read_csv($fn, ",");
    return $data[0];
}

sub read_nextclade_data{
    my $fn = shift;
    my @data = read_csv($fn, '\t');
    return $data[0];
}

sub read_csv {
    my $fn = shift;
    my $sep = shift;
    open (my $fh, $fn);
    chomp(my $header = <$fh>);
    $header =~ s/\r//;
    my @header = split /$sep/, $header;
    
    my @data;
    while(<$fh>) {
	chomp;
	s/\r//;
	my @a = split /$sep/;
	my %entry;
	for my $i (0..$#header) {
	    $entry{$header[$i]} = $a[$i];
	}
	push @data, \%entry;
    }
    return @data;
    
}

sub read_fasta_dumb {
    my $fn = shift;
    open(my $fh, $fn);
    my $header = <$fh>;
    chomp( my $seq = <$fh>);
    return $seq;
}

sub read_qual_dumb {
    my $fn = shift;
    open(my $fh, $fn);
    chomp( my $qual = <$fh>);
    return $qual;
}
    

sub get_ontarget_from_flagstat {
    my $fn = shift;
    open(my $fh, $fn);
    while(<$fh>) {
	chomp;
	s/\r//;
	if (/mapped \((.*?)\% :/) {
	    return $1;
	}
    }	
    return "N/A"
}

sub is_indel {
    my $v = shift;
    my ($pos,$bases) = split /_/, $v;
    my ($ref,$alt) = split '>', $bases;
    return 1 if length($ref) != length($alt);
    return 0;
}

sub read_depth_data {
    my $fn = shift;
    my $sample_id = shift;
    my $sample_oid = shift;
    my $variants = shift;

    my %positions;
    for my $v (@$variants) {
	my $pos = (split /_/, $v)[0];
	$positions{$pos}=1;

	# For indels, add the position after as well
	if( is_indel($v) ) {
	    $positions{$pos+1}=1;
	}
    }

    open (my $fh, $fn);
    <$fh>;
    my @depth;
    while (my $line=<$fh>) {
	chomp $line;
	my @a = split /\t/, $line;
	if( $positions{$a[1]+1} ) {
	    my %depth_obj = ('sample_oid'=>$sample_oid, 'sample_id'=>$sample_id, 'pos'=>$a[1]+1, 'dp'=>$a[2], 'A'=>$a[3], 'C'=>$a[4], 'G'=>$a[5], 'T'=>$a[6], 'DEL'=>$a[7], 'REFSKIP'=>$a[8]);
	    push @depth, \%depth_obj;
	}
    }
    return @depth;
}



sub get_depth_data_for_variants {
    my $dir = shift;
    my $variants_old = shift;
    my $variants_new = shift;
    my $samples_old = shift;
    my $samples_new = shift;

    my @all_variants = (@$variants_old, @$variants_new);

    my @all_depth_data;
    
    # Add depth data for new variants to all old samples
    for my $oid ( keys %{$samples_old} ) {
	my @dp = read_depth_data($samples_old->{$oid}->{depth_filename}, $samples_old->{$oid}->{sample_id}, $oid, $variants_new);
	push(@all_depth_data, @dp);
    }


    # Add depth data for all variants for the new samples
    for my $sample ( @$samples_new ) {
	my $sid = $sample->{sample_id};
	my $soid = $sample->{sample_oid};
	my @dp = read_depth_data($sample->{depth_filename}, $sid, $soid, \@all_variants);
	push(@all_depth_data, @dp);
    }

    return @all_depth_data;
}

sub get_inserted_sample_oids {
    my $insert_res = shift;
    my $samples = shift;
 
    my %oids;
    for my $inserted (@{$insert_res->{inserted}}) {
	my $sample_oid = $inserted->{'_id'}->value;
	my $idx = $inserted->{index};
	$oids{$samples->[$idx]->{sample_id}} = $sample_oid;
    }
    return %oids;
}

sub add_sample_oid {
    my $oids = shift;
    my $data = shift;

    for( my $i=0; $i<scalar(@$data); $i++ ) {
	my $sample_id = $data->[$i]->{sample_id};
	if($oids->{$sample_id}) {
	    $data->[$i]->{sample_oid} = $oids->{$sample_id};
	}
	else {
	    die "Something went wrong. No matching samples object ID för $sample_id\n";
	}
    }
}
