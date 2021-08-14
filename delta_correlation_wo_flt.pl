#!/usr/bin/perl
use Math::CDF;
open WT,"/home/lingshi/correlation/temp/wt.TPM.txt";
open IDH,"/home/lingshi/correlation/temp/mut.TPM.txt";
$line=<WT>;
chomp $line;
@info=split(/\t/,$line);
$n2=scalar(@info)-1;
$line=<IDH>;
chomp $line;
@info=split(/\t/,$line);
$n1=scalar(@info)-1;
open WTRS,"/home/lingshi/correlation/metadata/same_result_180kb_wt_TPM_wo_flt.txt";
open SAMED,">/home/lingshi/correlation/metadata/same_delta_180kb_TPM_wo_flt.txt";
$line=<WTRS>;
while($line=<WTRS>){
    chomp $line;
    @info=split(/\t/,$line);
    open IDHRS,"/home/lingshi/correlation/metadata/same_result_180kb_mut_TPM_wo_flt.txt";
    $line1=<IDHRS>;
    while($line1=<IDHRS>){
        chomp $line1;
        @info1=split(/\t/,$line1);
        if((($info[0] eq $info1[0] && $info[1] eq $info1[1]) || ($info[0] eq $info1[1] && $info[1] eq $info1[0])) && $info[5] < 1 && $info1[5] < 1){
            $d=$info1[5]-$info[5];
            $z2=0.5*log((1+$info[5])/(1-$info[5]));
            $z1=0.5*log((1+$info1[5])/(1-$info1[5]));
            $z=($z1-$z2)/sqrt(1/($n1-3)+1/($n2-3));
            $p= 2*(1-(Math::CDF::pnorm(abs($z))));
            print SAMED "$info[0]\t$info[1]\t$d\t$p\n";
        }
    }
    close IDHRS;
}
open WTRD,"/home/lingshi/correlation/metadata/dif_result_180kb_wt_TPM_wo_flt.txt";
#open IDHRD,"/home/lingshi/correlation/metadata/dif_result_180kb_mut_TPM_wo_flt.txt";
open DIFD,">/home/lingshi/correlation/metadata/dif_delta_180kb_TPM_wo_flt.txt";
$line=<WTRD>;
while($line=<WTRD>){
    chomp $line;
    @info=split(/\t/,$line);
    open IDHRD,"/home/lingshi/correlation/metadata/dif_result_180kb_mut_TPM_flt.txt";
    $line1=<IDHRD>;
    while($line1=<IDHRD>){
        chomp $line1;
        @info1=split(/\t/,$line1);
        if((($info[0] eq $info1[0] && $info[1] eq $info1[1]) || ($info[0] eq $info1[1] && $info[1] eq $info1[0])) && $info[5] < 1 && $info1[5] < 1){
            $d=$info1[5]-$info[5];
            $z2=0.5*log((1+$info[5])/(1-$info[5]));
            $z1=0.5*log((1+$info1[5])/(1-$info1[5]));
            $z=($z1-$z2)/sqrt(1/($n1-3)+1/($n2-3));
            $p= 2*(1-(Math::CDF::pnorm(abs($z))));
            print DIFD "$info[0]\t$info[1]\t$d\t$p\n";
        }
    }
    close IDHRD;
}
