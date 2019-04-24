#!/usr/bin/perl -w
# Modified perl script to process alignment data. 
#
use strict;
use warnings;
use Carp qw| cluck :DEFAULT |;
use Getopt::Long qw(:config gnu_getopt auto_version auto_help);
use Pod::Usage;
BEGIN { our $VERSION = "1.0"; }

our $opt;
if( !GetOptions( $opt = { debug => 0, }, 'debug!', 'man!', 'query|q=s', 'quiet!' ) ) { pod2usage(2); }

my $dbg = $opt->{debug};
if($opt->{man}){ pod2usage(-verbose => 2); }

if(!$opt->{query}) { pod2usage(  -message => "Error: required --query argument is missing!\n", -exitval => 1, -verbose => 1, -output  => \*STDERR ); }

my $query;
my $qstart;
my $seqsinfile = _read_fasta_file();
my @seqs;

# $seq = [ 'QUERY/1-126', 'ADKELKFLVVDDFSTMRRIVRNLLKELGFNNVEEAEDG...' ]
foreach my $seq (@$seqsinfile)
{
    if(!defined($query) && $seq->[0] =~ /$opt->{query}/o)
    {
        $qstart = ( defined($1) ? $1+0 : undef );
        unshift @seqs, ( $query = $seq ); # query may match multiple times, take the first only
    }
    else
    {
        push @seqs, $seq;
    }
}

my @querycols = ();
for( my $i = 0, my $s = $query->[1], my $l = length($query->[1]); $i < $l; ++$i )
{
    if(substr($s, $i, 1) =~ /[A-Z]/o){ push @querycols, $i; }
    #  if(substr($s, $i, 0) =~ /[>]/o){ push @querycols, $i; }
}

#if(defined($qstart)) { print STDOUT "# querystart=$qstart\n"; }
#{
#    my $qnogaps = $query->[1]; $qnogaps =~ tr/[A-Za-z]//cd;
#    print STDOUT "# query=", $qnogaps, "\n";
#}
#
#
foreach my $seq (@seqs)
{

    if(length($query->[1]) != length($seq->[1]))
    {
        warn("length of aligned sequence '".$seq->[0]."' (".length($seq->[1]).") does not equal to query length (".length($query->[1])."), skipping sequence\n");
        next;
    }
    # we just do default METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES=2 of calculate_evolutionary_constraints.m: skip entire sequence
    #if($seq->[1] =~ /[BJOXZ]/o){ if(!$opt->{quiet}){ warn("skipped '".$seq->[0]."' because it has [BJOXZ] in it\n"); } next; }

    print STDOUT ">.$seq->[0]\n";
    print STDOUT @{[split(//o, $seq->[1])]}[ @querycols ], "\n"; # neat
}

exit(0);

sub                 _read_fasta_file
{
    my $ret = [];
    my( $desc, $seq ) = ( undef, '' );
    my $line = <STDIN>;
    while($line)
    {
        if(substr($line,0,1) ne '>'){ $line = <STDIN>; next; }
        #
        if(defined($desc)){ push @$ret, [ $desc, $seq ]; $desc = undef; $seq = ''; }
        #
        chomp($line);
        $desc = substr($line, 1);
        $line = <STDIN>;
        while($line)
        {
            if(substr($line,0,1) ne '>'){ chomp($line); $seq .= $line; $line = <STDIN>; }
            else{ last; }
        }
    }
    if(defined($desc)){ push @$ret, [ $desc, $seq ]; }
    #
    return $ret;
}

=pod

=head1 NAME

a2m2aln - reformat A2M input to a simple alignment format

=head1 SYNOPSIS

a2m2aln [OPTIONS]

a2m2aln --query '^RASH_HUMAN/(\d+)' --quiet < INFILE > OUTFILE

a2m2aln --man --help --version

=head1 DESCRIPTION

a2m2aln formats L<A2M|http://compbio.soe.ucsc.edu/a2m-desc.html> input to a simple alignment format used by freecontact(1):

 * Optional header line: '# querystart=<position>'.
 * Optional header line: '# query=<SEQUENCE>'. Lowercase letters in <SEQUENCE> indicate insertions that were deleted from the alignment.
 * One aligned sequence per line.
 * The first sequence is the query.

All gaps and insertions - also query insertions - are removed from the alignment.  The 'query' header field helps reconstruct original query residue numbers.

=head1 OPTIONS

=over

=item -q, --query REGEXP

Query identifier, a regular expression, e.g. '^RASH_HUMAN\b' to match 'RASH_HUMAN/5-165'.  Required.

Use parentheses to match query start, e.g. '^RASH_HUMAN/(\d+)'. Matching the query start position is optional.

=item --debug

=item --nodebug

=item --help

=item --man

=item --quiet

Suppress printing of messages like 'sequences skipped for unusual residue letters'.

=item --version

=back

=head1 AUTHOR

Laszlo Kajan <lkajan@rostlab.org>

=head1 SEE ALSO

L<freecontact(1)>

=cut

#:vim:ts=4:et:ai:
