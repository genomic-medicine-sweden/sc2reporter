#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;


# Connect to database, and create handles for collections
my $client = MongoDB->connect();
my $SAMPLE = $client->ns("sarscov2_standalone.sample");


my $db_add_date = $ARGV[0];
my $ct_file = $ARGV[1];

open(CT, $ct_file) or die "cannot open file";
my $in_results = 0;
my %ct;
while(<CT>) {
    chomp;
    if( /^\[Results\]/ ) {
	<CT>;
	$in_results = 1;
	next;
    }
    if($in_results) {
	my @a = split /\t/;
	$ct{$a[0]} = $a[2];
    }
}
close CT;

my ($y, $m, $d) = split /-/, $db_add_date;

my $start_time = DateTime->new( year => $y, month => $m, day => $d );
my $end_time = $start_time->clone();
$end_time->add(days => 1);

foreach my $id ( keys %ct) {
    my $add_metadata_ref = $SAMPLE->update_one({'sample_id'=>$id, 'time_added'=>{'$gte'=>$start_time, '$lte'=>$end_time}}, {'$set'=>{'Ct'=>$ct{$id}}});
}

 
