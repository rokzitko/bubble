#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use Math::Trig qw(pi);

my $bubble = "$Bin/../bubble";
my $params = "$Bin/param.loop";
my $reference = "$Bin/cond.opt.geo-mma.dat";
my $resigma = "$Bin/resigma.dat";
my $imsigma = "$Bin/imsigma.dat";
my $output = "$Bin/cond.opt.geo.dat";
my $absolute_tolerance = 1e-12;
my $relative_tolerance = 3e-5;
my $quiet = 0;
my $write_output = 1;

if (@ARGV == 1 && $ARGV[0] eq '--quiet') {
    $quiet = 1;
} elsif (@ARGV == 1 && $ARGV[0] eq '--check') {
    $quiet = 1;
    $write_output = 0;
} elsif (@ARGV != 0) {
    die "Usage: $0 [--quiet|--check]\n";
}

-x $bubble or die "Executable not found: $bubble\n";

open my $param_file, '<', $params or die "Cannot open $params: $!\n";
my $temperature;
while (my $line = <$param_file>) {
    if ($line =~ /^\s*T\s*=\s*([^\s#]+)/) {
        $temperature = $1;
        last;
    }
}
close $param_file or die "Cannot close $params: $!\n";
defined $temperature or die "Parameter T not found in $params\n";

my $number = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;
$temperature =~ /^$number$/ or die "Invalid T in $params: $temperature\n";

open my $reference_file, '<', $reference or die "Cannot open $reference: $!\n";
my @mesh;
while (my $line = <$reference_file>) {
    next if $line =~ /^\s*$/;
    my ($omega, $expected) = split ' ', $line;
    defined $omega && defined $expected && $omega =~ /^$number$/ && $expected =~ /^$number$/
        or die "Malformed reference row: $line";
    push @mesh, [$omega, $expected + 0.0];
}
close $reference_file or die "Cannot close $reference: $!\n";
@mesh == 50 or die "Expected 50 optical frequencies in $reference, found " . scalar(@mesh) . "\n";

my @results;
my $passed = 0;
my $failed = 0;
my ($max_absolute_error, $max_relative_error) = (0.0, 0.0);
for my $index (0 .. $#mesh) {
    my ($omega, $expected) = @{$mesh[$index]};
    my @command = (
        $bubble,
        '-i', '2',
        '-k', '1',
        '-a', '1e-13',
        '-r', '1e-8',
        '-c', '20',
        '-s', '1e-16',
        '-O', $omega,
        '5', '2', '0', $temperature, '0', $resigma, $imsigma,
    );

    open my $pipe, '-|', @command or die "Cannot run $bubble: $!\n";
    my $raw = do { local $/; <$pipe> };
    close $pipe or die "bubble failed for Omega=$omega with exit status " . ($? >> 8) . "\n";
    defined $raw && $raw =~ /^\s*($number)\s*$/
        or die "Unexpected bubble output for Omega=$omega: " . (defined $raw ? $raw : '<empty>') . "\n";

    my $conductivity = 2*pi*($1 + 0.0);
    my $absolute_error = abs($conductivity - $expected);
    my $relative_error = $expected == 0.0 ? 0.0 : $absolute_error/abs($expected);
    my $allowed_error = $absolute_tolerance + $relative_tolerance*abs($expected);
    my $ok = $absolute_error <= $allowed_error;
    $ok ? ++$passed : ++$failed;
    $max_absolute_error = $absolute_error if $absolute_error > $max_absolute_error;
    $max_relative_error = $relative_error if $relative_error > $max_relative_error;

    push @results, [$omega, $conductivity];
    unless ($quiet) {
        printf "%2d/%d\t%s\t%.16g\t%s\trel=%.6g\n",
            $index + 1, scalar(@mesh), $omega, $conductivity, $ok ? 'OK' : 'FAILED', $relative_error;
    }
}

if ($write_output) {
    open my $output_file, '>', $output or die "Cannot open $output: $!\n";
    for my $row (@results) {
        printf {$output_file} "%s\t%.16g\n", $row->[0], $row->[1];
    }
    close $output_file or die "Cannot close $output: $!\n";
}

print "Wrote $output\n" if $write_output && !$quiet;
printf "BETHE OPTICAL OK=%d FAILED=%d max_abs=%.6g max_rel=%.6g\n",
    $passed, $failed, $max_absolute_error, $max_relative_error;
exit($failed == 0 ? 0 : 1);
