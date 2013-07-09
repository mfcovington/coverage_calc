package coverage_commander;
use Moose;
# use MooseX::UndefTolerant;
use Modern::Perl;
use Parallel::ForkManager;
use autodie;
# use Data::Printer;
use Capture::Tiny 'capture_stderr';
use CoverageDB::Main;
use POSIX;
use namespace::autoclean;

extends 'generic';

#TODO: check for presence of valid region!!!!
#TODO: require certain arguments to be defined
#TODO: generate log files??
#will it cause a problem if i look up a region that has no coverage?  will it return empty string, undef or 0?....looks like empty or undef
#TODO: Do I need to make defaults lazy and uncomment UndefTolerant?

sub samtools_cmd_gaps {
    my $self = shift;

    my $samtools_cmd = "samtools mpileup" . $self->_region . $self->bam . " | cut -f1-2,4 > " . $self->out_file . ".cov_gaps";
    return $samtools_cmd;
}

sub samtools_cmd_nogaps {
    my $self = shift;

    my $samtools_cmd = "samtools depth" . $self->_region . $self->bam . " > " . $self->out_file . ".cov_nogaps";
    return $samtools_cmd;
}

sub get_coverage {
    my $self = shift;

    $self->_validity_tests();
    $self->_make_dir();

    if ( $self->gap ) {
        say "  Running: " . $self->samtools_cmd_gaps() if $self->verbose();
        system( $self->samtools_cmd_gaps );
    }
    if ( $self->nogap ) {
        say "  Running: " . $self->samtools_cmd_nogaps() if $self->verbose();
        system( $self->samtools_cmd_nogaps );
    }
}

sub get_coverage_all {
    my $self = shift;

    if ( $self->gap ) {
        say "  Running: " . $self->samtools_cmd_gaps() if $self->verbose();
        system( $self->samtools_cmd_gaps );
    }
    if ( $self->nogap ) {
        say "  Running: " . $self->samtools_cmd_nogaps() if $self->verbose();
        system( $self->samtools_cmd_nogaps );
    }
}

around 'get_coverage_all' => sub {
    my $orig = shift;
    my $self = shift;

    $self->_validity_tests();

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager($self->threads);
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->_chromosome($chr);
        $self->out_file( $self->out_dir . "/coverage/" . $self->id . "." . $self->_chromosome . ".coverage" );
        $self->_make_dir();

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};

has 'cov_pos' => (
    is  => 'rw',
    isa => 'HashRef'
);

has 'cov_data' => (
    is  => 'ro',
    isa => 'ArrayRef'
);

has 'flank_dist' => (
    is      => 'rw',
    isa     => 'Int',
    default => 8,
);

has 'pos_start' => (
    is  => 'rw',
    isa => 'Int',
);

has 'pos_end' => (
    is  => 'rw',
    isa => 'Int',
);

has 'db' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 1,
);

has 'gap' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 1,
);

has 'nogap' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 1,
);

sub add_positions {
    my $self = shift;

    my $chr = $self->_chromosome;
    my $flank_dist = $self->flank_dist;

    my %cov_pos; # = $self->cov_pos;

# TODO: custom path
    open my $snps_fh, "<", "../genotyping/snp_master/polyDB.$chr.nr";
    <$snps_fh>;
    while (<$snps_fh>) {
        my $snp_pos = [ split /\t/ ]->[1];
        $cov_pos{$chr}{$snp_pos}                 = 1;
        $cov_pos{$chr}{ $snp_pos - $flank_dist } = 1;
        $cov_pos{$chr}{ $snp_pos + $flank_dist } = 1;
    }
    close $snps_fh;
    # print scalar keys $cov_pos{$chr}, "\n";
    $self->cov_pos( \%cov_pos );
}

sub get_coverage_db {
    my $self = shift;

    $self->add_positions;
    $self->populate_CoverageDB_by_chr;
}

around 'get_coverage_db' => sub {
    my $orig = shift;
    my $self = shift;

    $self->_validity_tests();

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager( floor $self->threads / 2 );
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->_chromosome($chr);

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};

sub populate_CoverageDB_by_chr {
    my $self = shift;

# TODO: custom path (and make empty db?)
    my $dbi = 'SQLite';
    my $db = 'db/coverage.db';
    my $schema = CoverageDB::Main->connect("dbi:$dbi:$db");

    my $chromosome  = $self->_chromosome;
    my $flank_dist  = $self->flank_dist;
    my $cov_pos_ref = $self->cov_pos;
    my $bam_file    = $self->bam;
    my $sample_id   = $self->id;

    say "  Getting coverage for $chromosome" if $self->verbose;
    my $count = 1;
    my @cov_data;

    system("samtools index $bam_file") if ! -e "$bam_file.bai";

    my $sam_gap_cmd = "samtools mpileup -r $chromosome $bam_file | cut -f1-2,4";
    my $sam_nogap_cmd = "samtools depth -r $chromosome $bam_file";

    my $gap_fh;
    capture_stderr {    # suppress mpileup output sent to stderr
        open $gap_fh,   "-|", $sam_gap_cmd;
    };
    open my $nogap_fh, "-|", $sam_nogap_cmd;
    while ( my $gap_line = <$gap_fh> ) {
        my $nogap_line = <$nogap_fh>;
        chomp( $gap_line, $nogap_line );
        my ( $chr, $pos, $gap_cov ) = split /\t/, $gap_line;
        my $nogap_cov = [ split /\t/, $nogap_line ]->[2];
        if ( exists $$cov_pos_ref{$chr}{$pos} ) {
            $count++;
            push @cov_data, [ $sample_id, $chr, $pos, $gap_cov, $nogap_cov ];
        }

        populate_and_reset( \$count, \@cov_data, \$schema ) if $count % 100000 == 0;
    }
    close $gap_fh;
    close $nogap_fh;

    populate_and_reset( \$count, \@cov_data, \$schema ) if scalar @cov_data;
}

sub populate_and_reset {
    my ( $count_ref, $cov_data_ref, $schema_ref ) = @_;
    $$count_ref = 1;
    $$schema_ref->populate(
        'Coverage',
        [
            [qw/sample_id chromosome position gap_cov nogap_cov/],
            @$cov_data_ref
        ]
    );
    @$cov_data_ref = ();
}

sub _region {
    my $self = shift;

    my $region;
    given ( $self->_chromosome ) {
        $region = " -r " . $self->_chromosome . ":" . $self->pos_start . "-" . $self->pos_end . " "
            when defined
            and defined $self->pos_start
            and defined $self->pos_end
            and $self->pos_start < $self->pos_end;
        $region = " -r " . $self->_chromosome . " " when defined;
        default { $region = " " }
    }

    return $region;
}

__PACKAGE__->meta->make_immutable;
