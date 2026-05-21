#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(mkpath);
use File::Basename;
use File::Spec;
use File::Find;

# 定义命令行参数
my ($outputdir, $inputdir, $threads, $max_parallel);

sub parse_command_line {
    GetOptions(
        "outputdir|o=s" => \$outputdir,          # 输出目录
        "inputdir|i=s"  => \$inputdir,           # 输入目录
        "threads=i"     => \$threads,            # 使用的线程数
        "maxparallel=i" => \$max_parallel,       # 最大并行任务数
    ) or die "无法解析命令行参数：$!";

    die "输出目录未指定\n" unless defined $outputdir;
    die "输入目录未指定\n" unless defined $inputdir;
}

sub set_default_values {
    $max_parallel  ||= 5;
    $threads       ||= 4;
}

parse_command_line();
set_default_values();

mkpath($outputdir) unless -d $outputdir;

# 递归查找 fastq 文件
my @samples;
find(sub {
    return unless -f $_;
    return unless /\.(fastq|fq)$/;
    push @samples, $File::Find::name;
}, $inputdir);

if (scalar @samples == 0) {
    die "未找到任何 .fastq 文件在输入目录: $inputdir\n";
}

my @active_processes;
my $qc_report_dir = "$outputdir/qc_reports";
mkpath($qc_report_dir) unless -d $qc_report_dir;

my %processed_samples;

foreach my $file (@samples) {
    my ($sample_id, $read_type);
    
    # 获取相对于inputdir的相对路径，例如 "homo/rawdata/GSE263588/SRR..."
    my $rel_path = File::Spec->abs2rel($file, $inputdir);
    my $rel_dir = dirname($rel_path);
    
    # 移除 'rawdata' 层级
    # 例如: "homo/rawdata/GSE263588" -> "homo/GSE263588"
    my $clean_rel_dir = $rel_dir;
    $clean_rel_dir =~ s{/rawdata/}{/}g;
    $clean_rel_dir =~ s{^rawdata/}{}g;
    $clean_rel_dir =~ s{/rawdata$}{}g;
    
    if ($file =~ /(.+)_1\.(fastq|fq)$/) {
        $sample_id = basename($1);
    } elsif ($file =~ /(.+)_2\.(fastq|fq)$/) {
        $sample_id = basename($1);
    } else {
        next;
    }
    
    # 唯一标识符
    my $unique_id = "$clean_rel_dir/$sample_id";
    next if $processed_samples{$unique_id};
    $processed_samples{$unique_id} = 1;

    print "开始质控样本: $unique_id\n";
    
    # 使用清理后的路径结构
    my $target_out_dir = "$outputdir/$clean_rel_dir";
    mkpath($target_out_dir) unless -d $target_out_dir;

    my $fq1 = "$inputdir/$rel_dir/${sample_id}_1.fastq";
    my $fq2 = "$inputdir/$rel_dir/${sample_id}_2.fastq";
    my $fq1_out = "$target_out_dir/${sample_id}_1.clean.fastq";
    my $fq2_out = "$target_out_dir/${sample_id}_2.clean.fastq";
    
    # 报告也按目录存放
    my $report_subdir = "$qc_report_dir/$clean_rel_dir";
    mkpath($report_subdir) unless -d $report_subdir;
    my $html_report = "$report_subdir/${sample_id}_fastp.html";
    my $json_report = "$report_subdir/${sample_id}_fastp.json";

    my $pid = fork();
    if (!defined $pid) {
        die "无法创建子进程: $!";
    } elsif ($pid == 0) { 
        my $cmd;
        if (-e $fq1 && -e $fq2) {
            $cmd = "fastp -i $fq1 -I $fq2 -o $fq1_out -O $fq2_out -h $html_report -j $json_report -w $threads";
        } elsif (-e $fq1) {
            $cmd = "fastp -i $fq1 -o $fq1_out -h $html_report -j $json_report -w $threads";
        } else {
            exit(1);
        }
        system($cmd) == 0 or die "fastp 执行失败\n";
        exit(0);
    } else {
        push @active_processes, $pid;
    }

    while (@active_processes >= $max_parallel) {
        my $child = waitpid(-1, 0);
        @active_processes = grep { $_ != $child } @active_processes;
    }
}

while (wait() != -1) { }

my $finished_flag = "$outputdir/QC_finished.txt";
open(my $FH, '>', $finished_flag) or die "无法创建标志文件: $!";
print $FH "All QC jobs finished.\n";
close($FH);
