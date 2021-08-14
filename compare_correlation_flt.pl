#!/usr/bin/perl
$prefix=$ARGV[0];
$pvalue=$ARGV[1];
$same_wt=$prefix.'_wt_TPM_flt.txt';
$same_mut=$prefix.'_mut_TPM_flt.txt';
$prefix=~s/same/dif/g;
$dif_wt=$prefix.'_wt_TPM_flt.txt';
$dif_mut=$prefix.'_mut_TPM_flt.txt';
open SAME_WT,"$same_wt";
open SAME_MUT,"$same_mut";
open DIF_WT,"$dif_wt";
open DIF_MUT,"$dif_mut";
$same_cpr=$same_wt;
$same_cpr=~s/result/compare/g;
$same_cpr=~s/_wt//g;
open SAME_CPR,">$same_cpr";
$dif_cpr=$same_cpr;
$dif_cpr=~s/same/dif/g;
open DIF_CPR,">$dif_cpr";
$line=<SAME_WT>;
while($line=<SAME_WT>){
    chomp $line;
    @info=split(/\t/,$line);
    $id=$info[0]."\t".$info[1];
    $same{$id}=$info[5];
    $samep{$id}=$info[6];
}
$line=<SAME_MUT>;
print SAME_CPR "geneid1\tgeneid2\tgene1\tgene2\tdistance\tcorrelation_change\tp_value\tTAD\n";
while($line=<SAME_MUT>){
    chomp $line;
    @info=split(/\t/,$line);
    $id1=$info[0]."\t".$info[1];
    $id2=$info[1]."\t".$info[0];
    if(exists $same{$id1}){
        $fold=$info[5]-$same{$id1}; 
        $q=$info[6]*$samep{$id1};
        if($q<$pvalue){
            print SAME_CPR "$id1\t$info[2]\t$info[3]\t$info[4]\t$fold\t$q\t$info[7]\n";
        }
    }
    elsif(exists $same{$id2}){
        $fold=$info[5]-$same{$id2};
        $q=$info[6]*$samep{$id2};
        if($q<$pvalue){
            print SAME_CPR "$id2\t$info[2]\t$info[3]\t$info[4]\t$fold\t$q\t$info[7]\n";
        }
    }
}
$line1=<DIF_WT>;
while($line1=<DIF_WT>){
    chomp $line1;
    @info1=split(/\t/,$line1);
    $id=$info1[0]."\t".$info1[1];
    $dif{$id}=$info1[5];
    $difp{$id}=$info1[6];
}
$line1=<DIF_MUT>;
print DIF_CPR "geneid1\tgeneid2\tgene1\tgene2\tdistance\tcorrelation_change\tp_value\tTAD1\tTAD2\n";
while($line1=<DIF_MUT>){
    chomp $line1;
    @info1=split(/\t/,$line1);
    $id1=$info1[0]."\t".$info1[1];
    $id2=$info1[1]."\t".$info1[0];
    if(exists $dif{$id1}){
        $fold=$info1[5]-$dif{$id1};
        $q=$info1[6]*$difp{$id1};
        if($q<$pvalue){
            print DIF_CPR "$id1\t$info1[2]\t$info1[3]\t$info1[4]\t$fold\t$q\t$info1[7]\t$info1[8]\n";
        }
    }
    elsif(exists $dif{$id2}){
        $fold=$info1[5]-$dif{$id2};
        $q=$info1[6]*$difp{$id2};
        if($q<$pvalue){
            print DIF_CPR "$id2\t$info1[2]\t$info1[3]\t$info1[4]\t$fold\t$q\t$info1[7]\t$info1[8]\n";
        }
    }
}
