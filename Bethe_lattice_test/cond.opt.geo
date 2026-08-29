#!/usr/bin/env perl

use strict;
use warnings;
use Errno qw(ERANGE);
use FindBin qw($Bin);
use IO::Select;
use IPC::Open3;
use Math::Trig qw(pi);
use POSIX qw(isfinite strtod WNOHANG);
use Symbol qw(gensym);
use Time::HiRes qw(time);

my $bubble = "$Bin/../bubble";
my $params = "$Bin/param.loop";
my $reference = "$Bin/cond.opt.geo-mma.dat";
my $resigma = "$Bin/resigma.dat";
my $imsigma = "$Bin/imsigma.dat";
my $output = "$Bin/cond.opt.geo.dat";
# The frozen Mathematica interpolant clips knots but not inter-knot overshoot.
my $absolute_tolerance = 3e-11;
my $relative_tolerance = 3e-5;
my $quiet = 0;
my $write_output = 1;
my $number = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;
my $command_timeout = $ENV{BUBBLE_TEST_TIMEOUT} // 120;
my $output_limit = 1024*1024;

if (@ARGV == 1 && $ARGV[0] eq '--quiet') {
    $quiet = 1;
} elsif (@ARGV == 1 && $ARGV[0] eq '--check') {
    $quiet = 1;
    $write_output = 0;
} elsif (@ARGV != 0) {
    die "Usage: $0 [--quiet|--check]\n";
}

-x $bubble or die "Executable not found: $bubble\n";

sub parse_finite_number {
    my ($token) = @_;
    return undef unless defined $token && $token =~ /\A$number\z/;
    local $! = 0;
    my ($value, $unparsed) = strtod($token);
    return undef if $! == ERANGE || $unparsed != 0 || !isfinite($value);
    return $value;
}

sub capture_command {
    my (@command) = @_;
    my $stdout = gensym;
    my $stderr = gensym;
    my $pid = open3(undef, $stdout, $stderr, @command);
    my %captured = (stdout => '', stderr => '');
    my %stream_name = (
        fileno($stdout) => 'stdout',
        fileno($stderr) => 'stderr',
    );
    my $selector = IO::Select->new($stdout, $stderr);
    my $deadline = time + $command_timeout;
    my $capture_error;
    while ($selector->count) {
        my $remaining = $deadline - time;
        if ($remaining <= 0) {
            $capture_error = "command timed out after $command_timeout seconds";
            last;
        }
        my @ready = $selector->can_read($remaining);
        if (!@ready) {
            $capture_error = "command timed out after $command_timeout seconds";
            last;
        }
        for my $handle (@ready) {
            my $fd = fileno($handle);
            my $bytes = sysread($handle, my $chunk, 8192);
            if (!defined $bytes) {
                $capture_error = "unable to read $stream_name{$fd}: $!";
                last;
            }
            if ($bytes == 0) {
                $selector->remove($handle);
                close $handle;
                next;
            }
            $captured{$stream_name{$fd}} .= $chunk;
            if (length($captured{stdout}) + length($captured{stderr}) >
                $output_limit) {
                $capture_error = "command output exceeded $output_limit bytes";
                last;
            }
        }
        last if defined $capture_error;
    }

    my $wait_status;
    unless (defined $capture_error) {
        while (1) {
            my $waited = waitpid($pid, WNOHANG);
            if ($waited == $pid) {
                $wait_status = $?;
                last;
            }
            if ($waited == -1) {
                $capture_error = "unable to wait for command: $!";
                last;
            }
            my $remaining = $deadline - time;
            if ($remaining <= 0) {
                $capture_error = "command timed out after $command_timeout seconds";
                last;
            }
            my $pause = $remaining < 0.05 ? $remaining : 0.05;
            select undef, undef, undef, $pause;
        }
    }
    if (defined $capture_error) {
        kill 'KILL', $pid;
        waitpid $pid, 0;
        die "Test harness error: $capture_error\n";
    }

    my $status = ($wait_status & 127)
        ? 128 + ($wait_status & 127) : $wait_status >> 8;
    return ($status, $captured{stdout}, $captured{stderr});
}

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
defined parse_finite_number($temperature)
    or die "Invalid T in $params: $temperature\n";

open my $reference_file, '<', $reference or die "Cannot open $reference: $!\n";
my @mesh;
while (my $line = <$reference_file>) {
    next if $line =~ /^\s*$/;
    my @fields = split ' ', $line;
    @fields == 2 or die "Malformed reference row: $line";
    my ($omega, $expected_token) = @fields;
    my $omega_value = parse_finite_number($omega);
    my $expected = parse_finite_number($expected_token);
    defined $omega_value && defined $expected
        or die "Malformed reference row: $line";
    push @mesh, [$omega, $expected];
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
        '-E',
        '-q',
        '-O', $omega,
        '5', '2', '0', $temperature, '0', $resigma, $imsigma,
    );

    my ($status, $raw, $stderr) = capture_command(@command);
    $status == 0
        or die "bubble failed for Omega=$omega with exit status $status: $stderr\n";
    $stderr eq ''
        or die "bubble emitted unexpected stderr for Omega=$omega: $stderr\n";
    defined $raw && $raw =~ /^\s*($number)\s*$/
        or die "Unexpected bubble output for Omega=$omega: " . (defined $raw ? $raw : '<empty>') . "\n";

    my $value = parse_finite_number($1);
    defined $value
        or die "Unrepresentable bubble output for Omega=$omega: $raw\n";
    my $conductivity = 2*pi*$value;
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
